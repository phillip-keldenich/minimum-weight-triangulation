import os, sys, random, re, json
import slurminade
import subprocess
import socket

__file__ = os.path.abspath(__file__)

slurminade.update_default_configuration(
    partition="alg",
    constraint="alggen05",
    mail_user="p.keldenich@tu-braunschweig.de",
    mail_type="FAIL",
    exclusive=True
)


uniform_normal_name_re = re.compile(r"^(uniform|normal)-([0-9]+)-([0-9]+)$")


def gurobi_license_file():
    hostname = socket.gethostname().split(".")[0]
    if not hostname.startswith("alg"):
        raise RuntimeError(f"HOSTNAME: '{hostname}' does not look as expected!")
    return os.path.join(os.environ["HOME"], ".gurobi", f"{hostname}-gurobi.lic")


def instance_name(filename):
    if filename.endswith(".bz2"):
        return instance_name(filename[:-4])
    if filename.endswith(".tsp") or filename.endswith(".txt"):
        return filename[:-4]
    if filename.endswith(".instance"):
        return filename[:-9]
    if filename.endswith(".instance.json"):
        return filename[:-14]
    return None


class ExperimentHandler:
    def __init__(self, num_repeats, only_inexact):
        self.num_repeats = num_repeats
        self.only_inexact = only_inexact
        self.instance_data_loc = os.path.abspath(os.path.dirname(__file__) + "/../../../instance_sets")
        self.experiment_wrapper_loc = os.path.abspath(os.path.dirname(__file__) + "/../../experiment_wrapper.py")
        self.executable_loc = os.path.abspath(os.path.dirname(__file__) + "/../../build/Release/src/solve_mwt")
        self.verifier_loc = os.path.abspath(os.path.dirname(__file__) + "/../../build/Release/src/verify_triangulation")
        if not os.path.isfile(self.experiment_wrapper_loc):
            raise RuntimeError(f"Experiment wrapper script missing (expected location: {self.experiment_wrapper_loc})")
        if not os.path.isfile(self.executable_loc):
            raise RuntimeError(f"Release-built executable missing (expected location: {self.executable_loc})")
        if not os.path.isfile(self.verifier_loc):
            raise RuntimeError(f"Release-built verifier missing (expected location: {self.verifier_loc})")
        self.instance_sets = [x for x in os.listdir(self.instance_data_loc) if os.path.isdir(os.path.join(self.instance_data_loc, x))]
        self.output_data_loc = os.path.join(os.path.dirname(__file__), "02_raw_output")
        if not os.path.isdir(self.output_data_loc):
            try:
                os.mkdir(self.output_data_loc)
            except FileExistsError:
                pass
        self.all_instances = []
        for instance_set in self.instance_sets:
            try:
                os.mkdir(os.path.join(self.output_data_loc, instance_set))
            except FileExistsError:
                pass
            instance_set_dir = os.path.join(self.instance_data_loc, instance_set)
            for f in os.listdir(instance_set_dir):
                name = instance_name(f)
                if name is None:
                    continue
                m = uniform_normal_name_re.match(name)
                if not m or int(m[2]) <= 50000000:
                    self.all_instances.append((name, instance_set, os.path.join(instance_set_dir, f)))

    def instance_output_prefix(self, instance_tuple):
        iname, iset, _ = instance_tuple
        return os.path.join(self.output_data_loc, iset, f"{iname}")

    def generate_commands_for(self, instance_tuple):
        out_prefix = self.instance_output_prefix(instance_tuple)
        modes = ("exact", "inexact")
        stdout_paths = [out_prefix + "_" + m + "_run{run}.stdout.log" for m in modes]
        stderr_paths = [out_prefix + "_" + m + "_run{run}.stderr.log" for m in modes]
        output_paths = [out_prefix + "_" + m + "_run{run}.json" for m in modes]
        mode_flags = {"exact": "--nonsimple-face-mode=exact",
                      "inexact": "--nonsimple-face-mode=inexact_with_gap"}
        result = []
        for run in range(1, self.num_repeats + 1):
            for i in range((1 if self.only_inexact else 0), 2):
                errpath = stderr_paths[i].format(run=run)
                outpath = stdout_paths[i].format(run=run)
                if os.path.isfile(errpath) or os.path.isfile(outpath) or os.path.isfile(errpath + ".bz2") or os.path.isfile(outpath + ".bz2"):
                    print(f"Skipping existing run '{outpath}'; clear the run to rerun it")
                    continue
                command = [
                    sys.executable, self.experiment_wrapper_loc,
                    "-o", outpath,
                    "-e", errpath,
                    "--tee", f"--time-limit={90 * 60}", "--print-end-message", 
                    "--"
                ]
                mode = modes[i]
                command += [self.executable_loc, mode_flags[mode], "--print-step-times", "--lp-verbose", "--monitor-resources",
                            "--verification=none", "-o", output_paths[i].format(run=run), "-i", instance_tuple[2]]
                result.append(command)
        return result


@slurminade.slurmify()
def run_command_group(command_group, verifier_loc):
    newenv = os.environ.copy()
    lfile = gurobi_license_file()
    newenv["GRB_LICENSE_FILE"] = lfile
    for command in command_group:
        # run command (blocks until termination)
        subprocess.run(command, capture_output=True, text=True, check=False, env=newenv)
        # run verification
        for i, p in enumerate(command):
            if p == '-o' and i + 1 < len(command):
                path = command[i + 1]
                if os.path.isfile(path) and path.endswith(".json"):
                    res = subprocess.run(['/usr/bin/time', verifier_loc, '-i', path, '--xfree-mode=sweep', '--print-steps'],
                                         capture_output=True, check=False)
                    verification_data = {"verification_stdout": str(res.stdout, encoding='utf-8'),
                                         "verification_stderr": str(res.stderr, encoding='utf-8'),
                                         "verification_returncode": res.returncode,
                                         "verification_successful": res.returncode == 0,
                                         "verification_input_path": path}
                    verdatapath = path[:-5] + "_verification.json"
                    with open(verdatapath, "w") as verfile:
                        json.dump(verification_data, verfile)
                    subprocess.run(["bzip2", "--best", verdatapath], capture_output=True, check=False)
        # post-processing: zip text outputs with bzip2
        for i, p in enumerate(command):
            if p == '-e' or p == '-o':
                if i + 1 < len(command):
                    path = command[i + 1]
                    if os.path.isfile(path):
                        subprocess.run(["bzip2", "--best", path], capture_output=True, check=False)


if __name__ == "__main__":
    NUM_REPEATS = 5
    ONLY_INEXACT = False
    slurminade.set_dispatch_limit(10000)
    e = ExperimentHandler(NUM_REPEATS, ONLY_INEXACT)
    all_commands = []
    for instance in e.all_instances:
        all_commands += e.generate_commands_for(instance)
    random.Random().shuffle(all_commands)

    answer = input(f"Found {len(all_commands)} commands to run -- continue? [y/n]")
    if answer.lower() != 'y' and answer.lower() != 'z':
        print("Aborting...\n")
        sys.exit(0)

    # reduce the number of elements pushed to slurm 
    # by making internal groups of 10 commands each
    job_bundle_size = 400
    command_group_size = 1
    command_groups = [[]]
    for command in all_commands:
        this_group = command_groups[-1]
        if len(this_group) >= command_group_size:
            command_groups.append([])
            this_group = command_groups[-1]
        this_group.append(command)

    with slurminade.JobBundling(max_size=job_bundle_size):
        for command_group in command_groups:
            run_command_group.distribute(command_group, e.verifier_loc)

