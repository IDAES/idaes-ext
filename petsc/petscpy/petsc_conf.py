import os
import json

def get_conf():
     pth = os.path.join(os.path.dirname(__file__), "petsc_conf.json")
     with open(pth, "r") as fp:
         cfg = tuple(json.load(fp))
     return cfg
