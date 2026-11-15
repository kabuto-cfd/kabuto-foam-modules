# OpenFOAM-14 in the Claude workspace

The Claude workspace (Ubuntu 24.04, amd64) reaches GitHub but not openfoam.org, so it cannot
install the `openfoam14` package itself. The workflow `openfoam14-runtime.yml` repackages that
package unchanged as a release asset of this repository (`openfoam14-runtime-<version>`), with
the matching source pack, `DEPENDS.txt` and `BUILD-INFO.txt`.

In the workspace `gh release download` fails (it uses GraphQL, which is blocked); the REST API
works, and the 120 MB archive comes down in seconds:

```bash
repo=kabuto-cfd/kabuto-foam-modules
rel=$(gh api repos/$repo/releases --jq \
  '[.[] | select(.tag_name | startswith("openfoam14-runtime-"))][0].tag_name')
id=$(gh api repos/$repo/releases/tags/$rel --jq \
  '.assets[] | select(.name | endswith("-ubuntu24.04-amd64.tar.xz")) | .id')
gh api -H "Accept: application/octet-stream" repos/$repo/releases/assets/$id > /tmp/openfoam14.tar.xz
sudo tar -C /opt -xJf /tmp/openfoam14.tar.xz
sudo apt-get install -y -qq --no-install-recommends g++ make flex binutils-dev \
    libreadline-dev libxt-dev zlib1g libopenmpi-dev openmpi-bin
export USER=${USER:-root}       # else the user directories become ~/OpenFOAM/-14
set +e; source /opt/openfoam14/etc/bashrc; set -e
foamRun -help > /dev/null && echo "OpenFOAM-14 ready"
```

The package brings its own Scotch (`/opt/openfoam14/thirdparty`) and the full source, so
`wmake` builds user libraries against it.

To refresh the release after an openfoam.org patch release, start the workflow by hand:

```bash
gh api -X POST repos/kabuto-cfd/kabuto-foam-modules/actions/workflows/openfoam14-runtime.yml/dispatches -f ref=main
```

(Pushing a tag from the workspace is refused by its git proxy; a tag pushed from elsewhere,
`openfoam14-runtime-<anything>`, also starts it.) The release is named after the package
version; the workflow fails if that release already exists.

If no release can be downloaded, build from source instead (2 to 4 hours on 2 cores; a
background command is stopped after 2 hours, rerun `./Allwmake -j 2 -q` to continue):

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
