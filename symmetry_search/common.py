from pathlib import Path
import json
import os
from sage.all import QQ, CyclotomicField, matrix

ROOT = Path(__file__).resolve().parent.parent

class LimitExceeded(RuntimeError):
    pass

def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(data, indent=2) + '\n')
    os.replace(tmp, path)

def field_from(conductor):
    return QQ if int(conductor) <= 2 else CyclotomicField(int(conductor), names='z')

def encode(x):
    if x in QQ:
        return [str(QQ(x))]
    return [str(c) for c in x.list()]

def decode(x, K):
    if K is QQ:
        if len(x) > 1:
            raise ValueError('Non-rational coefficient in rational model')
        return QQ(x[0]) if x else QQ(0)
    return K([QQ(c) for c in x])

def encode_matrix(M):
    return [[encode(x) for x in row] for row in M.rows()]

def decode_matrix(data,K):
    return matrix(K,[[decode(x,K) for x in row] for row in data])

def conductor(K):
    return 1 if K is QQ else int(K.zeta_order())
