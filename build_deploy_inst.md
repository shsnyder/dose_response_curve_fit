### Build (assuming the pyproject.toml file is in place)
- `python3 -m build`

### Deploy
- First time:
    - `pip install <path to wheel>`
- Update:
    - `pip install <path to wheel> --upgrade`

### Sample pyproject.toml contents
~~~text
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "skewnormal_curvefit"
version = "0.9"
authors = [
  { name="Scott Snyder ", email="sh_snyder@yahoo.com" },
]
description = "Curve fit of dose response data using Skew Normal distribution"
readme = "README.md"
requires-python = ">=3.7"
classifiers = [
    "Programming Language :: Python :: 3",
    "Operating System :: OS Independent",
]
~~~