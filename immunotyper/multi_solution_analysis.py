import os
from .models import ModelMeta
from gurobipy import Model as GurobiSolver
from .common import log, resource_path
from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Union, Any
import yaml
import json

class PrefixConsistencyError(Exception):
    """Custom exception for errors during prefix consistency calculation,
       typically due to missing prerequisite data like band solutions.
    """
    pass


class MultiSolutionAnalysis(object):

    def __init__(self, model: ModelMeta, gene_type: str, allele_db):
        """
        Initialize the MultiSolutionAnalysis object for finding multiple solutions
        in bands relative to the optimal objective function value.
        
        Args:
            model: A model object with a Gurobi solver
        
        Raises:
            ValueError: If the model does not use GurobiSolver
        """
        if not isinstance(model.model, GurobiSolver):
            raise ValueError("MultiSolutionAnalysis only supports GurobiSolver models")
        
        self.model = model
        self.gene_type = gene_type
        self.allele_db = allele_db
        self.gurobi_model = model.model
        self.allele_calls = model.get_allele_calls()
        self.functional_allele_calls = model.get_allele_calls(functional_only=True)
        self.objective_value = model.get_value(model.objective)
                
        self.optimal_solution = self._create_solution_info(0, self.objective_value)

        self.near_optimal_solutions, self.differing_alt_optimal = self.extract_optimal_solutions()
        print(f"Found {len(self.near_optimal_solutions)} near-optimal solutions")
        self.optimal_solutions = [s for s in self.near_optimal_solutions if s['is_optimal']]


    def get_multi_optimal_different_alleles(self, optimal_solutions=None) -> List[Dict[str, Any]]:
        """
        Get a list of optimal solutions that have different alleles than the original optimal solution.
        
        Returns:
            List of solution dictionaries containing solution information
        """
        if not optimal_solutions:
            optimal_solutions = self.optimal_solutions

        if not hasattr(self, 'optimal_solutions'):
            log.warning("No near_optimal_solutions found. Call extract_optimal_solutions() first.")
            return []
        
        # Get the optimal solution alleles
        optimal_alleles = set(self.optimal_solutions['allele_vars'].items())

        # Check optimal_solutions are actually optimal
        if not all([s['objective_value'] == self.objective_value for s in optimal_solutions]):
            raise ValueError("All optimal_solutions must have the same objective value as the original optimal solution")

        # Filter to get only optimal solutions with different alleles
        different_allele_solutions = []
        for solution in optimal_solutions:
            solution_alleles = set(solution['allele_vars'].items())
            if solution_alleles != optimal_alleles:
                different_allele_solutions.append(solution)
        
        log.info(f"Found {len(different_allele_solutions)} optimal solutions with different alleles")
        return different_allele_solutions

    def write_multiple_solutions_alleles(self, output_dir: str, output_prefix: str, optimal_solutions=None, only_differing=True) -> None:
        """
        Write each truly optimal solution to separate files in the specified directory if it has different allele calls from otpimal.
        
        For each optimal solution, two files are created:
        1. A file containing all allele calls (repeats each allele based on its copy number)
        2. A file containing only functional allele calls
        
        Args:
            output_dir: Directory where solution files will be written
            output_prefix: Prefix for output filenames
            
        Returns:
            None
            
        Note:
            This method requires that you've already called extract_optimal_solutions()
            and stored the results in self.near_optimal_solutions
            
        Files created:
            - {output_prefix}_optimal_{i}_allele_calls.txt: All allele calls for solution i
            - {output_prefix}_optimal_{i}_functional_allele_calls.txt: Functional allele calls for solution i
        """
        if not optimal_solutions:
            if only_differing:
                optimal = self.differing_alt_optimal
            else:
                optimal = self.optimal_solutions
        
        # Filter to get only truly optimal solutions
        optimal = self.differing_alt_optimal
        if not optimal:
            log.warning("No alternative optimal solutions, not writing any files")
            return
        
        log.info(f"Writing {len(optimal)} optimal solutions to {output_dir}")
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Process each optimal solution
        for i, solution in enumerate(optimal):
            allele_calls = []
            functional_allele_calls = []
            
            # Generate expanded allele lists (repeating each allele based on copy number)
            for allele, copy_number in solution['allele_vars'].items():
                for j in range(copy_number):
                    allele_calls.append(allele)
                    if self.allele_db[allele].is_functional:
                        functional_allele_calls.append(allele)
            
            # Write all allele calls
            allele_calls_path = f"{output_dir}/{output_prefix}-alt_optimal_{i}_allele_calls.txt"
            with open(allele_calls_path, "w") as f:
                f.write("\n".join(allele_calls))
            
            # Write functional allele calls
            functional_path = f"{output_dir}/{output_prefix}-alt_optimal_{i}_functional_allele_calls.txt"
            with open(functional_path, "w") as f:
                f.write("\n".join(functional_allele_calls))
            
            log.info(f"Wrote solution {i} with {len(allele_calls)} total alleles "
                    f"({len(functional_allele_calls)} functional) to {allele_calls_path}")
        
        if not optimal:
            log.warning("No truly optimal alternative solutions found within the pool.")
        else:
            log.info(f"Finished writing {len(optimal)} optimal solutions")

    def _create_solution_info(
        self, 
        solution_number: int, 
        objective_value: float,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Create a standardized solution information dictionary.
        
        Args:
            solution_number: Index number for the solution
            objective_value: Objective function value for this solution
            **kwargs: Additional fields to include in the solution info dict
        
        Returns:
            Dictionary containing standardized solution information
        """
        # Extract allele information and read errors from current model state
        allele_vars = self.extract_allele_information(self.model)
        read_error, read_error_norm, read_assignments = self.get_read_assignment_error()
        
        # Calculate percent gap from optimal
        percent_from_optimal = ((objective_value - self.objective_value) / self.objective_value) * 100
        
        # Create solution info with standard fields
        solution_info = {
            'solution_number': solution_number,
            'objective_value': objective_value,
            'percent_from_optimal': percent_from_optimal,
            'allele_vars': allele_vars,
            'read_errors': {'absolute': read_error, 'normalized': read_error_norm},
            'read_assignments': read_assignments
        }
        
        # Add any additional fields from kwargs
        solution_info.update(kwargs)
        
        return solution_info
        
    def get_read_assignment_error(self, model: Optional[ModelMeta] = None) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, List[Tuple[str, float]]]]:
        """
        Calculate read assignment errors for alleles in the model.
        
        Args:
            model: Model to analyze (defaults to self.model if None)
            
        Returns:
            Tuple containing:
                - allele_read_error: Dictionary mapping allele IDs to total read error
                - allele_read_error_norm: Dictionary mapping allele IDs to normalized read error
                - allele_read_assignment: Dictionary mapping allele IDs to list of (read_id, distance) tuples
        """
        if not model:
            model = self.model

        candidates = model._candidates
        alleles = model.get_allele_calls(functional_only=True)
        candidates_dict = dict([(c.id, c) for c in candidates])
        allele_read_assignment = defaultdict(list)
        allele_read_error = {} 
        allele_read_error_norm = {} # scaled by number of copies
        for a in alleles:
            # Get the candidate object for this allelez
            cand = candidates_dict[a]
            if cand.num_copies == 0:
                continue
            sum_error = 0
            for read in cand.reads:
                if model.get_value(cand.assignment[read.id]) == 1: # read was assigned to candidate
                    allele_read_assignment[cand.id].append((read.id, read.get_distance(cand)))
                    sum_error += read.get_distance(cand)
            allele_read_error[cand.id] = sum_error
            allele_read_error_norm[cand.id] = sum_error/cand.num_copies
        
        return allele_read_error, allele_read_error_norm, allele_read_assignment

    def extract_allele_information(self, model: ModelMeta) -> Dict[str, int]:
        """
        Extract allele copy number information from the current state of the model.
        
        Args:
            model: Model to extract allele information from
            
        Returns:
            Dictionary mapping allele IDs to their copy numbers
        """
        allele_vars = {}
        for c in model.candidates:
            if c.num_copies > 0:
                allele_vars[c.id] = c.num_copies
        return allele_vars

    def solve_multi_solution_bands(
        self, 
        num_bands: int = 5, 
        max_gap_percent: float = 10.0, 
        time_limit: int = 4000
    ) -> List[Dict[str, Any]]:
        """
        Solve for multiple solutions in bands relative to the optimal objective function value.
        Uses the optimal solution stored during initialization as band 0.
        
        Args:
            num_bands: Number of bands to divide the gap into
            max_gap_percent: Maximum gap percentage from optimal solution (e.g., 10 for 10%)
            time_limit: Time limit in seconds for each optimization
        
        Returns:
            List of solution dictionaries containing solution information
        """
        # Set time limit for optimization
        self.gurobi_model.params.timeLimit = time_limit
        
        # Get the optimal value
        optimal_value = self.objective_value
        log.info(f"Optimal solution value: {optimal_value}")
        
        # Initialize solutions list with the optimal solution
        all_solutions = []
        
        # Get a reference to the objective expression
        obj_expr = self.model.objective
        
        # Iterate through each band
        for band in range(num_bands):
            # Calculate band boundaries as percentages of optimal
            lower_pct = band * (max_gap_percent / num_bands) / 100
            upper_pct = (band +1) * (max_gap_percent / num_bands) / 100
            
            # Calculate absolute bounds for this band
            lower_bound = optimal_value * (1 + lower_pct)
            upper_bound = optimal_value * (1 + upper_pct)
            
            log.info(f"Band {band+1}: Finding solution with objective value between {lower_pct:.2f} and {upper_pct:.2f} difference from optimal")
            
            # Reset to clear solution information but keep parameters
            try:
                self.gurobi_model.reset(0)
            except:
                log.warning("Could not perform partial reset, continuing anyway")
            
            # Use standard optimization mode
            self.gurobi_model.setParam('PoolSearchMode', 0)
            
            # Add constraints for this band
            lower_bound_constr = self.gurobi_model.addConstr(obj_expr >= lower_bound, name=f"band_lower_{band}")
            upper_bound_constr = self.gurobi_model.addConstr(obj_expr <= upper_bound, name=f"band_upper_{band}")
            
            # Find solution in this band
            log.info(f"Optimizing for band {band+1}...")
            self.gurobi_model.optimize()
            
            # Check if we found a solution
            if self.gurobi_model.Status == 2:  # Optimal solution found
                obj_val = self.gurobi_model.objVal
                
                # Create solution info dictionary using helper function
                solution_info = self._create_solution_info(
                    solution_number=len(all_solutions),
                    objective_value=obj_val,
                    band=band
                )
                
                all_solutions.append(solution_info)
                
                log.info(f"Band {band} Solution: Objective Value = {obj_val:.2f} " +
                            f"({solution_info['percent_from_optimal']:.2f}% from optimal)")
            else:
                log.warning(f"No solution found for Band {band}. Status code: {self.gurobi_model.Status}")
            
            # Remove band constraints
            self.gurobi_model.remove([lower_bound_constr, upper_bound_constr])
        
        log.info(f"Total solutions found: {len(all_solutions)} across {num_bands} bands")
        self.multi_band_solutions = all_solutions

    def get_solutions_by_band(
        self,
        solutions: List[Dict[str, Any]],
        band: int
    ) -> List[Dict[str, Any]]:
        """
        Filter solutions to only those from a specific band.
        
        Args:
            solutions: List of solution dictionaries
            band: Band number to filter by (0 for optimal, 1-num_bands for banded solutions)
            
        Returns:
            List of solution dictionaries from the specified band
        """
        return [s for s in solutions if s['band'] == band]
        
    def extract_optimal_solutions(
        self,
        tolerance: float = 0.05  # Allow solutions within 0.01% of optimal
    ) -> List[Dict[str, Any]]:
        """
        Extract multiple optimal or near-optimal solutions from the Gurobi solution pool.
        
        Args:
            tolerance: Percentage tolerance from optimal to consider a solution "near-optimal" (default: 0.01%)
            
        Returns:
            List of solution dictionaries containing solution information
        """
        self.gurobi_model.Params.SolutionNumber = 0

        # Get the number of solutions in the pool
        sol_count = self.gurobi_model.SolCount
        log.debug(f"Found {sol_count} solutions in the Gurobi solution pool")
        
        if sol_count == 0:
            return []
        
        # Use the optimal solution we already stored
        best_obj_val = self.optimal_solution['objective_value']
        
        # Maximum acceptable objective value based on tolerance
        max_obj_val = best_obj_val * (1 + tolerance/100)
        
        # Initialize storage for solutions with the optimal solution
        pool_solutions = []
        
        # Iterate through all solutions in the pool
        for i in range(1, sol_count):
            self.gurobi_model.Params.SolutionNumber = i
            obj_val = self.gurobi_model.PoolObjVal
                        
            # Calculate gap percentage from optimal
            percent_gap = ((obj_val - best_obj_val) / best_obj_val) * 100
            
            # Determine if this solution is truly optimal (exact match)
            is_optimal = obj_val == best_obj_val

            # Only include solutions within the tolerance
            if obj_val <= max_obj_val:
                if is_optimal:
                    log.info(f"Solution {i}: Objective Value = {obj_val:.4f} (Optimal)")
                else:
                    log.info(f"Solution {i}: Objective Value = {obj_val:.4f} ({percent_gap:.4f}% from optimal)")
                
                # Update model variables to reflect this solution
                for c in self.model.candidates:
                    if hasattr(c, 'var') and c.var is not None:
                        var_name = c.var.VarName
                        c.num_copies = int(round(self.gurobi_model.getVarByName(var_name).Xn))
                
                # Create solution info dictionary using helper function
                solution_info = self._create_solution_info(
                    solution_number=len(pool_solutions),
                    objective_value=obj_val,
                    pool_index=i,
                    is_optimal=is_optimal
                )
                
                pool_solutions.append(solution_info)
            else:
                log.info(f"Solution {i}: Objective Value = {obj_val:.4f} ({percent_gap:.4f}% from optimal) - Exceeds tolerance")
        
        log.info(f"Extracted {len(pool_solutions)} solutions within {tolerance}% of optimal")
        alt_optimal = [s for s in pool_solutions if s['is_optimal']]
        differing_alt_optimal = []
        if len(alt_optimal) > 0:
            log.info(f"Found {len(alt_optimal)} alternative optimal solutions")
            differing_alt_optimal = self.get_multi_optimal_different_alleles(pool_solutions)
            log.info(f"Found {len(differing_alt_optimal)} alternative optimal solutions with different alleles")
        else:
            log.warning("No alternative optimal solutions found")

        self.gurobi_model.Params.SolutionNumber = 0

        return pool_solutions, differing_alt_optimal
    
    def write_solution_info(self, output_file: str, solutions: List[Dict[str, Any]]) -> None:
        """
        Write solution information to a JSON file.
        
        Args:
            output_file: Path to output file 
            solutions: List of solution dictionaries containing solution information
            
        Returns:
            None
        """
        
        with open(output_file, "w") as f:
            json.dump(solutions, f, indent=4)
        
        log.info(f"Wrote solution information to {output_file}")
    
    def write_all_solutions(self, output_dir: str, prefix: str) -> None:
        """
        Write all solutions to a YAML file.
        
        Args:
            output_file: Path to output file
            
        Returns:
            None
        """
        solution_dir = os.path.join(output_dir, f'{prefix}-multi_solution_info')
        os.makedirs(solution_dir, exist_ok=True)
        
        output_file = os.path.join(solution_dir, "optimal_solution.json")
        self.write_solution_info(output_file, self.optimal_solution)

        output_file = os.path.join(solution_dir, "multi_band_solutions.json")
        self.write_solution_info(output_file, self.multi_band_solutions)

        output_file = os.path.join(solution_dir, "near_optimal_solutions.json")
        self.write_solution_info(output_file, self.near_optimal_solutions)

    def calculate_prefix_consistency_for_optimal_alleles(self) -> Dict[str, Dict[str, Any]]:
        """
        Calculates the prefix consistency value for each allele call from the optimal solution.

        The prefix consistency score for an allele is defined as the number of consecutive
        bands immediately following the optimal solution where the allele's copy number
        remains identical to its copy number in the optimal solution.

        Raises:
            PrefixConsistencyError: If prerequisite data (optimal solution or multi-band solutions)
                                    is missing or invalid.

        Returns:
            Dict[str, Dict[str, Any]]: A dictionary where keys are allele IDs from the
            optimal solution. Values are dictionaries containing:
                'in_optimal': True
                'optimal_copies': Copy number in the optimal solution.
                'is_functional': Boolean indicating functionality.
                'prefix_consistency': The calculated prefix consistency score.
                'objective_diff_at_consistency_break': Percent objective difference from optimal
                                                       for the last band that maintained consistency.
                                                       None if prefix_consistency is 0.
                'band0_copies': Copy number in the optimal solution ("Consistency Band 0").
                'band1_copies' through 'band5_copies': Copy numbers in "Consistency Bands 1-5".
                                                        (0 if band not solved or allele absent).
        """
        # Step A: Prerequisites and Initialization
        if not hasattr(self, 'optimal_solution') or \
           not self.optimal_solution or \
           'allele_vars' not in self.optimal_solution:
            raise PrefixConsistencyError("Optimal solution data is not available or is invalid.")

        if not hasattr(self, 'multi_band_solutions') or self.multi_band_solutions is None:
            # Allow self.multi_band_solutions to be an empty list
            # If it's None, it means solve_multi_solution_bands was likely not called.
            raise PrefixConsistencyError(
                "Multi-band solutions attribute (`self.multi_band_solutions`) is None. "
                "Ensure solve_multi_solution_bands() has been run."
            )

        num_inferred_bands = len(self.multi_band_solutions)
        # Using a global 'log' or 'self.model.log' if available for warnings
        # For this standalone function, I'll assume a 'log' object exists (e.g., from 'from .common import log')
        # If using a class-specific or instance-specific logger:
        # logger = getattr(self.model, 'log', logging.getLogger(__name__)) # Example
        if num_inferred_bands == 0:
            # This situation is handled by the logic below; prefix_consistency will be 0.
            # log.warning( # Replace with actual logging mechanism
            print(
                "Warning: No solutions found in `multi_band_solutions`. "
                "Prefix consistency scores will be 0 for all alleles "
                "(indicating no consistency beyond the optimal solution itself)."
            )


        # Step B: Prepare Band Data Structure
        bands_data: Dict[int, Dict[str, Any]] = {}

        # "Consistency Band 0" (Data from Optimal Solution)
        bands_data[0] = {
            'alleles': self.optimal_solution['allele_vars'],
            'obj_diff': 0.0
        }

        # "Consistency Band 1" through "Consistency Band num_inferred_bands"
        for idx in range(num_inferred_bands):
            # Assuming self.multi_band_solutions[idx] is the solution for the idx-th band interval
            # (i.e., solution['band'] in the original creation was idx)
            current_band_solution = self.multi_band_solutions[idx]
            consistency_band_index = idx + 1
            bands_data[consistency_band_index] = {
                'alleles': current_band_solution.get('allele_vars', {}),
                'obj_diff': current_band_solution.get('percent_from_optimal', float('inf'))
            }

        # Step C & D: Calculate Metrics for Alleles in Optimal Solution
        prefix_consistency_results: Dict[str, Dict[str, Any]] = {}

        for allele_id, optimal_copy_number in self.optimal_solution['allele_vars'].items():
            if optimal_copy_number == 0:  # Should not happen if allele_vars is sparse
                continue

            metrics: Dict[str, Any] = {}
            metrics['in_optimal'] = True
            metrics['optimal_copies'] = optimal_copy_number
            
            if hasattr(self, 'allele_db') and allele_id in self.allele_db and \
               hasattr(self.allele_db[allele_id], 'is_functional'):
                metrics['is_functional'] = self.allele_db[allele_id].is_functional
            else:
                # Fallback or error if allele_db structure is not as expected
                # log.warning(f"Could not determine functionality for allele {allele_id} from self.allele_db")
                print(f"Warning: Could not determine functionality for allele {allele_id} from self.allele_db")
                metrics['is_functional'] = False # Default if lookup fails

            metrics['band0_copies'] = optimal_copy_number

            # Calculate prefix_consistency score
            current_prefix_consistency_score = 0
            if optimal_copy_number > 0:
                for j in range(1, num_inferred_bands + 1): # Iterate through Consistency Bands 1, 2, ...
                    cn_band_j = bands_data.get(j, {}).get('alleles', {}).get(allele_id, 0)
                    if cn_band_j == optimal_copy_number:
                        current_prefix_consistency_score += 1
                    else:
                        break  # Consistency with optimal_copy_number broken
            metrics['prefix_consistency'] = current_prefix_consistency_score

            # Determine objective_diff_at_consistency_break
            if metrics['prefix_consistency'] > 0:
                # The objective difference of the last band (j = metrics['prefix_consistency'])
                # that maintained consistency with the optimal solution.
                metrics['objective_diff_at_consistency_break'] = \
                    bands_data.get(metrics['prefix_consistency'], {}).get('obj_diff', None)
            else:
                metrics['objective_diff_at_consistency_break'] = None

            # Record copy numbers for a fixed set of bands (e.g., up to 5 beyond optimal)
            # This ensures a consistent structure in the output.
            max_bands_to_report = 5 # Report copy numbers for Consistency Bands 1 through 5
            for k in range(1, max_bands_to_report + 1):
                if k <= num_inferred_bands:
                    metrics[f'band{k}_copies'] = \
                        bands_data.get(k, {}).get('alleles', {}).get(allele_id, 0)
                else:
                    # This reporting band k is beyond the bands actually solved for/inferred
                    metrics[f'band{k}_copies'] = 0
            
            prefix_consistency_results[allele_id] = metrics

        return prefix_consistency_results


    def format_prefix_consistency_tsv(
        self,
        output_file_path: Optional[str] = None
    ) -> Optional[str]:
        """
        Formats the prefix consistency results for optimal alleles into a TSV string,
        integrating confidence metrics from an external JSON file.

        Args:
            confidence_metrics_file_path (str): The direct file system path to the
                                                confidence metrics JSON file.
            output_file_path (Optional[str]): If provided, the TSV output will be
                                              written to this file path. If None,
                                              the TSV string is returned.

        Returns:
            Optional[str]: A string formatted as a Tab-Separated Value (TSV) table
                           if output_file_path is None, otherwise None.

        Raises:
            FileNotFoundError: If the confidence_metrics_file_path does not exist.
            json.JSONDecodeError: If the confidence metrics file is not valid JSON.
            PrefixConsistencyError: If calculate_prefix_consistency_for_optimal_alleles()
                                    raises it due to its own prerequisites not being met.
            AttributeError: If self.gene_type is not set on the instance.
        """
        # 1. Resolve Path and Load External Confidence Metrics
        
        json_path = resource_path("confidence_metric/gene_prefix_consistency_metrics.json")
        
        try:
            with open(json_path, 'r') as f:
                confidence_data = json.load(f)
        except FileNotFoundError:
            log.error(
                f"Confidence metrics file not found at path: {json_path}"
            )
            return None
        except json.JSONDecodeError as e:
            log.error(
                f"Error decoding JSON from {json_path}: {e.msg} (Line: {e.lineno}, Col: {e.colno})"
            )
            return None

        # 2. Get Allele Prefix Consistency Data
        # This call might raise PrefixConsistencyError, which will propagate as per design.
        try:
            allele_consistency_data = self.calculate_prefix_consistency_for_optimal_alleles()
        except PrefixConsistencyError as e: # Explicitly catch to re-raise if needed, or handle
            log.error(f"Failed to calculate prefix consistency: {e}")
            raise # Re-raise if this error should halt execution further up.
                  # Or return None if this specific error should also be handled by returning None.
                  # For now, re-raising as it's a different category of error.


        # 3. Prepare for TSV Output (Store as list of dicts for sorting)
        processed_data_for_tsv: List[Dict[str, Any]] = []

        # 4. Iterate and Process Each Allele
        for allele_id, metrics_dict in allele_consistency_data.items():
            # a. Extract Basic Info
            functional_status_str = "Functional" if metrics_dict.get('is_functional', False) else "Non-functional"
            consistency_value = metrics_dict.get('prefix_consistency', 0)

            # b. Get gene_type (from instance) and base_gene_name
            if not hasattr(self, 'gene_type') or not self.gene_type:
                log.error("Instance attribute 'self.gene_type' is not set. Please ensure it's initialized.")
                return None # Or handle per allele: for now, fail the whole process
            current_gene_type = self.gene_type
            base_gene_name = allele_id.split('*')[0]

            # c. Initialize Derived Metrics to Default "Data missing"
            prob_tp_str = "Data missing"
            passes_optimal_str = "Data missing"
            optimal_thresh_accuracy_str = "Data missing"

            # d. Fetch Data from confidence_data
            gene_data_from_json = confidence_data.get("gene_types", {})\
                                               .get(current_gene_type, {})\
                                               .get("genes", {})\
                                               .get(base_gene_name, None)

            if gene_data_from_json:
                threshold_metrics_for_gene = gene_data_from_json.get("prefix_consistency_threshold_metrics", {})

                prob_tp_data_entry = threshold_metrics_for_gene.get(str(consistency_value), {})
                prob_tp_value = prob_tp_data_entry.get("precision", None)
                if prob_tp_value is not None:
                    prob_tp_str = f"{prob_tp_value:.2f}"

                optimal_threshold_val_from_json = gene_data_from_json.get("optimal_threshold", None)
                if optimal_threshold_val_from_json is not None:
                    passes = consistency_value >= optimal_threshold_val_from_json
                    passes_optimal_str = "True" if passes else "False"

                    optimal_thresh_key_str = str(optimal_threshold_val_from_json)
                    optimal_accuracy_data_entry = threshold_metrics_for_gene.get(optimal_thresh_key_str, {})
                    optimal_gene_accuracy_value = optimal_accuracy_data_entry.get("precision", None)
                    if optimal_gene_accuracy_value is not None:
                        optimal_thresh_accuracy_str = f"{optimal_gene_accuracy_value:.2f}"
            
            processed_data_for_tsv.append({
                "Allele": allele_id,
                "Functional status": functional_status_str,
                "Consistency value": consistency_value,
                "Probability of TP": prob_tp_str,
                "Passes Optimal Threshold": passes_optimal_str,
                "Optimal Threshold Gene Accuracy": optimal_thresh_accuracy_str
            })

        sorted_data_for_tsv = sorted(processed_data_for_tsv, key=lambda x: x["Allele"])

        header_columns = ["Allele", "Functional status", "Consistency value",
                          "Probability of TP", "Passes Optimal Threshold",
                          "Optimal Threshold Gene Accuracy"]
        tsv_rows: List[str] = ["\t".join(header_columns)]

        for row_dict in sorted_data_for_tsv:
            row_values = [
                str(row_dict["Allele"]),
                str(row_dict["Functional status"]),
                str(row_dict["Consistency value"]),
                str(row_dict["Probability of TP"]),
                str(row_dict["Passes Optimal Threshold"]),
                str(row_dict["Optimal Threshold Gene Accuracy"])
            ]
            tsv_rows.append("\t".join(row_values))

        tsv_string = "\n".join(tsv_rows)

        if output_file_path:
            try:
                with open(output_file_path, 'w') as f_out:
                    f_out.write(tsv_string)
                return None
            except IOError as e:
                log.error(f"Error writing TSV to file {output_file_path}: {e}")
                return None # Return None if file writing fails
        else:
            return tsv_string
