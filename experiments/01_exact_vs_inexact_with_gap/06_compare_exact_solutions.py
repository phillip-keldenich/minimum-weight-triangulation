import json, bz2, subprocess, os, re
from concurrent.futures import ProcessPoolExecutor


compare_program = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", "..", "build", "Release", "src", "compare_triangulations"))


def run_comparison_must_be_equal(runfile1, runfile2):
    result = subprocess.run([compare_program, "-i", runfile1, "-j", runfile2], text=True, check=False, capture_output=True)
    if result.returncode != 0:
        print("Error: comparison between", runfile1, runfile2)
        if result.stdout: print(result.stdout)
        if result.stderr: print(result.stderr)
        return False
    if "Interval difference: [0.0, 0.0]" not in result.stdout:
        print("Error: comparison not equal between", runfile1, runfile2)
        if result.stdout: print(result.stdout)
        if result.stderr: print(result.stderr)
        return False
    print(f"No error between '{runfile1}' and '{runfile2}'!")
    return True


interval_diff_line_re = re.compile(r"^\s*Interval\s+difference:\s+\[([0-9.eE+-]+),\s*([0-9.eE+-]+)\]\s*$")


def run_comparison_get_difference(better_file, worse_file):
    run_result = subprocess.run([compare_program, "-i", better_file, "-j", worse_file], text=True, check=False, capture_output=True)
    if run_result.returncode != 0:
        return {"error": True, "stdout": run_result.stdout, "stderr": run_result.stderr,
                "exit_code": run_result.returncode}
    result = {"error": False, "stdout": run_result.stdout, "stderr": run_result.stderr,
              "exit_code": run_result.returncode, "difference_lb": None, "difference_ub": None}
    for line in run_result.stdout.splitlines():
        if not line: continue
        m = interval_diff_line_re.match(line)
        if m is None: continue
        result["difference_lb"] = float(m[1])
        result["difference_ub"] = float(m[2])
        break
    if result["difference_lb"] is None:
        raise RuntimeError("Could not find interval difference in output of comparison program!")
    return result


def get_output_file(instance_name, entry, run_id, run_tag="exact"):
    basedir = os.path.dirname(__file__)
    instance_dir = entry["instance_dir"]
    return os.path.realpath(os.path.join(basedir, instance_dir, f"{instance_name}_{run_tag}_run{run_id}.json.bz2"))


if __name__ == "__main__":
    exact_solves_path = os.path.realpath(os.path.join(os.path.dirname(__file__), "05a_successful_exact_solves.json"))
    with open(exact_solves_path, "r") as f:
        successful_exact_runs = json.load(f)
    with ProcessPoolExecutor(4) as pool:
        for instance, entry in successful_exact_runs.items():
            runs = entry["runs"]
            if len(runs) == 1: continue
            baseline_run = runs[0]
            outfile1 = get_output_file(instance, entry, baseline_run)
            if not os.path.isfile(outfile1):
                raise RuntimeError(f"Run file '{outfile1}' does not exist!")
            for run_id in runs[1:]:
                outfile2 = get_output_file(instance, entry, run_id)
                if outfile1 == outfile2:
                    raise RuntimeError("Duplicate file!")
                if not os.path.isfile(outfile2):
                    raise RuntimeError("Output file", outfile2, "does not exist!")
                pool.submit(run_comparison_must_be_equal, outfile1, outfile2)

