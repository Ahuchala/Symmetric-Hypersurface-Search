"""Killable, reusable Sage processes. No forks of a live Sage session."""
import json
import os
import select
import subprocess
import sys
import tempfile
import time
from .common import ROOT


class PersistentWorker:
    def __init__(self,model,prime,*,algorithm='auto',method='jacobian',sparse=True):
        self.config=dict(model=str(model),prime=prime,algorithm=algorithm,method=method,sparse=sparse)
        self.process=None;self.stderr=None;self.buffer=b'';self.starts=0

    def _line(self,seconds):
        deadline=time.monotonic()+seconds
        while b'\n' not in self.buffer:
            remaining=deadline-time.monotonic()
            if remaining<=0 or not select.select([self.process.stdout],[],[],remaining)[0]:
                raise TimeoutError('Sage worker timed out')
            data=os.read(self.process.stdout.fileno(),65536)
            if not data:raise RuntimeError('Sage worker exited without a result')
            self.buffer+=data
        line,self.buffer=self.buffer.split(b'\n',1)
        return json.loads(line)

    def _start(self,seconds):
        self.stderr=tempfile.TemporaryFile()
        self.process=subprocess.Popen([sys.executable,'-m','symmetry_search.chart_worker'],cwd=ROOT,
            stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.stderr)
        self.starts+=1
        self.process.stdin.write((json.dumps(self.config)+'\n').encode());self.process.stdin.flush()
        if self._line(seconds)!={'ready':True}:raise RuntimeError('Invalid worker handshake')

    def check(self,pivot,seconds):
        start=time.perf_counter();stage='startup'
        try:
            if self.process is None:self._start(seconds)
            stage='chart'
            self.process.stdin.write((json.dumps(list(pivot))+'\n').encode());self.process.stdin.flush()
            result=self._line(seconds)
        except TimeoutError:
            result={'status':'timeout','seconds':seconds,'stage':stage,'pivot':list(pivot)}
            self.close()
        except Exception as exc:
            if self.process is not None and self.process.poll() is None:self.process.kill();self.process.wait()
            error=str(exc)
            if self.stderr is not None:
                self.stderr.seek(0);error+='\n'+self.stderr.read()[-3000:].decode(errors='replace')
            result={'status':'error','error':error,'stage':stage,'pivot':list(pivot)}
            self.close()
        result['worker_wall_seconds']=time.perf_counter()-start
        return result

    def close(self):
        if self.process is not None:
            if self.process.poll() is None:self.process.kill()
            self.process.wait()
            try:self.process.stdin.close()
            except BrokenPipeError:pass
            self.process.stdout.close()
            self.process=None
        if self.stderr is not None:self.stderr.close();self.stderr=None
        self.buffer=b''

    def __enter__(self):return self
    def __exit__(self,*exc):self.close()
