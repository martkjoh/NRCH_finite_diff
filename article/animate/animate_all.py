import subprocess
from multiprocessing import Pool


folder = "article/animate/"
def make_anim(i):
    subprocess.run("python3"+" "+folder+"animate_"+i+".py", shell=True)

vids = ["2", "5", "sol"]

with Pool(10) as pool:
    pool.map(make_anim, vids)
