# Example GRASS GIS Addon

This Repository contains a GRASS GIS addon including all meta files
like github workflows, linter configs, test setup etc.
When starting a new vector related GRASS GIS addon, it is recommended
to use this repository as a starting point.

See also common [GRASS GIS addon methods](https://github.com/mundialis/grass-gis-helpers) which can be reused when writing your own addon.

See also README about [How to create a GRASS GIS addon](https://github.com/mundialis/grass-gis-helpers/How-to-create-a-GRASS-GIS-addon.md) including
best practises, structure, how to name it and more sources which might be usefull.

### Use pre-commit

It is highly recommended to install and use [pre-commit](https://pre-commit.com)
before submitting any new or modified code or any other content. The pre-commit
Git hooks set checks validity and executes automated formatting for
a range of file formats, including Python. Pre-commit installs
all necessary tools in a virtual environment upon first use.

If you never used pre-commit before, you must start by installing it on your
system. You only do it once:

```bash
python -m pip install pre-commit
```

Pre-commit must then be activated in the code repository. Change the directory
to the root folder and use the `install` command:

```bash
cd <actinia-core_source_dir>

# once per repo
pre-commit install
```

Pre-commit will then be automatically triggered by the `git commit` command. If
it finds any problem it will abort the commit and try to solve it automatically.
In that case review the changes and run again `git add` and
`git commit`.

It is also possible to run pre-commit manually, e.g:

```bash
pre-commit run linting --all-files
```

Or to target a specific set of files:

```bash
pre-commit run --files src/*
```

The pre-commit hooks set is defined in
[.pre-commit-config.yaml](.pre-commit-config.yaml).

It is possible to temporally disable the pre-commit hooks in the repo, e.g. while
working on older branches:

```bash
pre-commit uninstall
```

And to reactivate pre-commit again:

```bash
git switch main
pre-commit install
```
