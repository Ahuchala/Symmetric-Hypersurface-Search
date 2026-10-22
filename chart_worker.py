"""Private persistent worker: one model/field per process, JSON-line protocol."""
import json
import sys
from .models import read_model
from .chart_engine import ChartChecker


def main():
    config=json.loads(sys.stdin.readline())
    checker=ChartChecker(read_model(config.pop('model')),**config)
    print(json.dumps({'ready':True}),flush=True)
    for line in sys.stdin:
        pivot=json.loads(line)
        try:result=checker.check(tuple(pivot))
        except Exception as exc:result={'status':'error','error':str(exc),'pivot':pivot}
        print(json.dumps(result),flush=True)

if __name__=='__main__':main()
