import sys, os, bz2, random, subprocess, json, re
from data_handling import list_runs_in, instance_type_directories, classify_run, load_instance_data


runname_re = re.compile("^([a-zA-Z0-9-_]+)_(i?n?exact)_run([0-9]+)")


def identify_failed_runs():
    instance_type_dirs = instance_type_directories()
    failed_runs = []
    for instance_type_dir in instance_type_dirs:
        runs_in = list_runs_in(instance_type_dir)
        for run in runs_in:
            classification = classify_run(instance_type_dir, run)
            if classification == "COMPLETED": continue
            failed_runs.append((instance_type_dir, run))
    return failed_runs

def identify_no_lp_needed():
    instance_type_dirs = instance_type_directories()
    no_lp_runs = []
    for instance_type_dir in instance_type_dirs:
        runs_in = list_runs_in(instance_type_dir)
        for run in runs_in:
            datafile = os.path.join(instance_type_dir, f"{run}.json.bz2")
            if not os.path.isfile(datafile): continue
            with bz2.open(datafile, "rt") as output:
                data = output.read(32768)
            if '"face_triangulation_stats": null' in data:
                no_lp_runs.append((instance_type_dir, run))
                print(instance_type_dir, run)
    return no_lp_runs

def identify_lp_runs(no_lp_runs):
    instance_type_dirs = instance_type_directories()
    lp_runs = []
    for instance_type_dir in instance_type_dirs:
        runs_in = list_runs_in(instance_type_dir)
        for run in runs_in:
            datafile = os.path.join(instance_type_dir, f"{run}.json.bz2")
            if not os.path.isfile(datafile): continue
            if run in no_lp_runs: continue
            lp_runs.append((instance_type_dir, run))
            print(instance_type_dir, run)
    return lp_runs

def _experiment_stat_stub_json(data_stream):
    data = data_stream.read(16384)
    begin = '"solution_info": {'
    index = data.find(begin)
    if index == -1: raise RuntimeError("Unexpected begin of file...")
    index += 17
    nest_depth = 1
    end = index + 1
    try:
        while nest_depth > 0:
            if data[end] == '{':
                nest_depth += 1
            elif data[end] == '}':
                nest_depth -= 1
            end += 1
        return json.loads(data[index:end])
    except IndexError:
        data += data_stream.read()
        jsondata = json.loads(data)
        return jsondata["solution_info"]
    
    
def extract_times_and_stats():
    result_data_list = []
    instance_data = load_instance_data()
    for instance_type_dir in instance_type_directories():
        runs_in = list_runs_in(instance_type_dir)
        for run in runs_in:
            result_data_list.append({})
            result_data = result_data_list[-1]
            m = runname_re.match(run)
            if not m:
                raise RuntimeError("Invalid run name " + run + "!")
            print(run)
            result_data["instance_name"] = m[1]
            result_data["instance_dir"] = instance_type_dir
            result_data["run_type"] = m[2]
            result_data["run_index"] = m[3]
            output_file = os.path.join(instance_type_dir, f"{run}.json.bz2")
            if not os.path.isfile(output_file):
                result_data["run_succeeded"] = False
                result_data["run_failure_type"] = classify_run(instance_type_dir, run)
            else:
                result_data["run_succeeded"] = True
                result_data["run_failure_type"] = None
                with bz2.open(output_file, "rt") as f:
                    result_data["run_stats"] = _experiment_stat_stub_json(f)
                    if "num_points" not in result_data["run_stats"]:
                        result_data["run_stats"]["num_points"] = instance_data[m[1]]["num_points"]
    return result_data_list


if __name__ == "__main__":
    brief_data = extract_times_and_stats()
    with open("brief_run_data.json", "w") as f:
        json.dump(brief_data, f, indent=2)

