import subprocess, json, os, re
from concurrent.futures import ThreadPoolExecutor

instance_set_info = json.load(open("../../../instance_sets/instance_info.json", "r"))
exe_path = "../../build/Release/src/solve_mwt"
ver_path = "../../build/Release/src/verify_triangulation"
large_uniform_re = re.compile("uniform-([0-9]+)-([0-9]+).txt.bz2")

if not os.path.isfile(exe_path):
    raise RuntimeError("Could not find solver executable!")


RUN_REPEAT = 5

for f in os.listdir("../../../instance_sets/random"):
    ifile = os.path.join("../../../instance_sets/random", f)
    m = large_uniform_re.match(f)
    if not m: continue
    s = int(m[1])
    if s <= 50_000_000: continue
    for mode in ("exact", "inexact"):
        for rep in range(1, RUN_REPEAT+1):
            run_name = f"uniform-{s}-{m[2]}_{mode}_run{rep}"
            outfile = os.path.join(".", "02_raw_output", "random", f"{run_name}.json")
            out = os.path.join(".", "02_raw_output", "random", f"{run_name}.stdout.log")
            err = os.path.join(".", "02_raw_output", "random", f"{run_name}.stderr.log")
            if os.path.isfile(f"{outfile}.bz2"): continue
            print("Starting run:", run_name)
            modeflag = ("exact" if mode == "exact" else "inexact_with_gap")
            proc1 = subprocess.run(["/usr/bin/time", "-f", 'Max RSS: %M KB', exe_path, "-i", ifile,
                                    "--print-step-times", f"--nonsimple-face-mode={modeflag}", "--verification=none",
                                    "-o", outfile, "--monitor-resources", "--lp-verbose"], check=True, capture_output=True, text=True)
            print("Verifying result:", outfile)
            proc2 = subprocess.run([ver_path, "-i", outfile], check=True, capture_output=True, text=True)
            with open(out, "w") as of:
                of.write(proc1.stdout)
            with open(err, "w") as of:
                of.write(proc1.stderr)
            print("Compressing output...")
            with ThreadPoolExecutor(3) as pool:
                pool.submit((lambda: subprocess.run(["bzip2", "--best", outfile], check=True)))
                pool.submit((lambda: subprocess.run(["bzip2", "--best", out], check=True)))
                pool.submit((lambda: subprocess.run(["bzip2", "--best", err], check=True)))
