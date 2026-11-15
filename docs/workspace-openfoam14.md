# OpenFOAM-14 in the Claude workspace

The Claude workspace (Ubuntu 24.04, amd64) reaches GitHub but not openfoam.org, so it cannot
install the `openfoam14` package. The workflow `openfoam14-runtime.yml` repackages that package
unchanged as a release asset of this repository, with the matching source pack:

```bash
tag=$(gh release list --repo kabuto-cfd/kabuto-foam-modules --json tagName \
      --jq '[.[] | select(.tagName | startswith("openfoam14-runtime-"))][0].tagName')
mkdir -p /tmp/of14 && gh release download "$tag" --repo kabuto-cfd/kabuto-foam-modules \
      --pattern 'openfoam14-*-ubuntu24.04-amd64.tar.xz' --dir /tmp/of14
sudo tar -C /opt -xJf /tmp/of14/openfoam14-*.tar.xz
sudo apt-get install -y -qq $(grep -oP 'Depends: \K[^<>]+' DEPENDS.txt)   # runtime libraries
set +e; source /opt/openfoam14/etc/bashrc; set -e
foamRun -help > /dev/null && echo "OpenFOAM-14 ready"
```

To refresh it after an openfoam.org patch release, push a tag `openfoam14-runtime-<anything>`
(or run the workflow by hand); the release is named after the package version.

If the release cannot be downloaded, build from source instead (2 to 4 hours on 2 cores):

```bash
sudo apt-get install -y build-essential flex bison zlib1g-dev libopenmpi-dev openmpi-bin \
    libscotch-dev libptscotch-dev libreadline-dev libncurses-dev
mkdir -p ~/OpenFOAM && cd ~/OpenFOAM
git clone --depth 1 --branch <patch tag, e.g. 20261008> https://github.com/OpenFOAM/OpenFOAM-14
cd OpenFOAM-14
sed -i 's/^export SCOTCH_TYPE=ThirdParty/export SCOTCH_TYPE=system/;
        s/^export ZOLTAN_TYPE=ThirdParty/export ZOLTAN_TYPE=none/' etc/bashrc
set +e; source etc/bashrc; set -e
./Allwmake -j 2 -q > ../build14.log 2>&1
```

A background command in the workspace is stopped after 2 hours: rerun `./Allwmake -j 2 -q`,
it continues where it stopped.
