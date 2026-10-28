"""Run the full analysis in order. Steps share one namespace, so run them through this file."""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.join(HERE, 'analysis')
sys.path.insert(0, ANALYSIS)

STEPS = sorted(f for f in os.listdir(ANALYSIS) if f[:2].isdigit() and f.endswith('.py'))

ns = {'__name__': '__main__'}
exec(compile(open(os.path.join(ANALYSIS, 'common.py')).read(), 'common.py', 'exec'), ns)
for step in STEPS:
    t0 = time.time()
    print(f"\n[{step}]")
    path = os.path.join(ANALYSIS, step)
    exec(compile(open(path).read(), path, 'exec'), ns)
    print(f"[{step} done in {time.time() - t0:.0f}s]")
