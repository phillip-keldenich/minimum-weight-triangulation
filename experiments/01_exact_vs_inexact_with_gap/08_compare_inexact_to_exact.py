import json, bz2, subprocess, os, re
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Lock


compare_exact = __import__("06_compare_exact_solutions")
run_comparison_get_difference = compare_exact.run_comparison_get_difference
lock = Lock()

def comparison_to_dbfile(exact_run_file, inexact_run_file):
    result = run_comparison_get_difference(exact_run_file, inexact_run_file)
    result.update({"exact_run": exact_run_file, "inexact_run": inexact_run_file})
    new_line = json.dumps(result)
    with lock:
        print(new_line)


if __name__ == "__main__":
    exact_solves_path = os.path.realpath(os.path.join(os.path.dirname(__file__), "05a_successful_exact_solves.json"))
    with open(exact_solves_path, "r") as f:
        successful_exact_runs = json.load(f)
    with ProcessPoolExecutor(6) as pool:
        for instance, entry in successful_exact_runs.items():
            exact_run_file = os.path.join(entry["instance_dir"], f"{instance}_exact_run{entry['runs'][0]}.json.bz2")
            for inexact_run_id in range(1, 6):
                inexact_run_file = os.path.join(entry["instance_dir"], f"{instance}_inexact_run{inexact_run_id}.json.bz2")
                if not os.path.isfile(inexact_run_file):
                    continue
                pool.submit(comparison_to_dbfile, exact_run_file, inexact_run_file)
        
