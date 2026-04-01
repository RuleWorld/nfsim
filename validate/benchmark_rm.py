import time
import subprocess
import os

def run_bench(args):
    start = time.time()
    subprocess.check_call(args, stdout=subprocess.DEVNULL)
    return time.time() - start

print("Benchmarking v09.xml (High-rejection ring closure model)...")

nfsim_bin = os.path.join("build", "NFsim")
model = os.path.join("validate", "basicModels", "v09.xml")

time_standard = run_bench([nfsim_bin, "-xml", model, "-sim", "10", "-seed", "123"])
time_rm = run_bench([nfsim_bin, "-xml", model, "-sim", "10", "-seed", "123", "-rulemonkey"])

print(f"Standard NFsim Time: {time_standard:.3f}s")
print(f"RuleMonkey Time: {time_rm:.3f}s")

if time_rm < time_standard:
    print("RuleMonkey is faster for this high-rejection model! \n")
else:
    print("RuleMonkey is slower, which may happen depending on the model's exact parameters.")
