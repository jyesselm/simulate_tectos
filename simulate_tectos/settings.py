import os

file_path = os.path.realpath(__file__)
spl = file_path.split("/")
BASE_DIR = "/".join(spl[:-1])
RESOURCE_DIR = BASE_DIR + "/resources/"
