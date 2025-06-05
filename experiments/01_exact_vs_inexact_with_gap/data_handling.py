import os, re, bz2, json


runfile_re = re.compile(r"([a-zA-Z0-9_-]+)_(i?n?exact)_run([0-9]+)\.(stdout\.log|stderr\.log|json)(\.bz2)?")


def instance_type_directories():
    res = []
    local_dir = "./02_raw_output"
    for x in os.listdir(local_dir):
        p = os.path.join(local_dir, x)
        if os.path.isdir(p):
            res.append(p)
    return res


def list_runs_in(directory):
    run_names = set()
    for f in os.listdir(directory):
        m = runfile_re.match(f)
        if not m:
            if "verification" not in f:
                print("No match:", f)
            continue
        run_name = m[1] + "_" + m[2] + "_run" + m[3]
        run_names.add(run_name)
    return run_names


def classify_run(directory, run_name):
    main_name = os.path.join(directory, run_name)
    if os.path.isfile(main_name + ".json.bz2"):
        return "COMPLETED"
    outlog = main_name + ".stdout.log.bz2"
    errlog = main_name + ".stderr.log.bz2"
    if not os.path.isfile(errlog) or\
       not os.path.isfile(outlog):
        return "TASK_DIED"
    with bz2.open(errlog, "rt") as errf:
        error_output = errf.read()
    if "Gurobi seems to be unavailable" in error_output:
        return "GUROBI_UNAVAILABLE"
    if error_output.strip() == "" or "ommand terminated by signal 9" in error_output:
        return "TIMEOUT"
    if "Failed to read instance file" in error_output:
        return "INVALID_INSTANCE"
    if "Inexact solution is not close to integral" in error_output:
        return "EXACT_HAS_NO_IP"
    if "std::bad_alloc" in error_output:
        return "BAD_ALLOC"
    return "UNKNOWN"


def load_instance_data():
    path = os.path.join(os.path.dirname(__file__), "..", "..", "..", "instance_sets", "instance_info.json")
    with open(path, "r") as f:
        return json.load(f)


def classify_runs_in(directory):
    run_names = list_runs_in(directory)
    result = {}
    for run_name in run_names:
        result.setdefault(classify_run(directory, run_name), []).append(run_name)
    return result


def clean_task(directory, run_name):
    main_name = os.path.join(directory, run_name)
    for ext in (".json", ".stdout.log", ".stderr.log"):
        ext_name = main_name + ext
        if os.path.isfile(ext_name):
            os.remove(ext_name)
        if os.path.isfile(ext_name + ".bz2"):
            os.remove(ext_name + ".bz2")

