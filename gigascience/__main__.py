"""
Entry point to run GigaScience analysis module
"""

from subprocess import run

# run(["python", "-m", "gigascience.harmony"])
run(["python", "-m", "gigascience.regress"])
run(["python", "-m", "gigascience.cellxgene.py"])
