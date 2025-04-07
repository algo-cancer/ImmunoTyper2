import os
from .models import ModelMeta
from gurobipy import Model as GurobiSolver
from .common import log
from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Union, Any
import yaml
import json


class MultiSolutionAnalysis(object):

    def __init__(self, model: ModelMeta):
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

