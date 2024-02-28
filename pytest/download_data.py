"""
Download SatKit test vectors
to a local directory
"""

import json
from pathlib import Path
import requests
from download_from_json import download_from_json

if __name__ == "__main__":
    baseurl = "https://storage.googleapis.com/astrokit-astro-data"
    basedir = "astro-data"
    fileurl = baseurl + "/files.json"
    headers = {"Accept": "application/json"}

    Path(basedir).mkdir(exist_ok=True, parents=True)

    data = requests.get(fileurl, headers=headers).json()
    download_from_json(data, basedir, baseurl)
