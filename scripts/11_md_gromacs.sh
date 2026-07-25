#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────────
# Estágio 11 — dinâmica molecular do construto (GROMACS).
#
# Por que rodar MD depois do AlphaFold: o modelo predito é uma foto estática, e o
# construto é uma proteína quimérica cheia de linkers flexíveis — exatamente o tipo
# de molécula em que "a estrutura" não é uma estrutura só. O que a MD responde é se
# o novelo se mantém compacto ou se desenrola, e se os epitopos ficam expostos.
#
# Três réplicas com sementes diferentes, não uma só. Uma trajetória única de 100 ns
# não distingue relaxamento real de um artefato da velocidade inicial sorteada — e
# reportar RMSD de réplica única é o erro mais comum nesse tipo de artigo.
#
# Análises: RMSD (estabilidade global), RMSF por resíduo (quais blocos são móveis —
# linkers DEVEM ser móveis), raio de giro (compactação) e ligações de hidrogênio
# intramoleculares (integridade do enovelamento).
#
# Uso:
#     bash scripts/11_md_gromacs.sh 8
#     THREADS=16 bash scripts/11_md_gromacs.sh
# ─────────────────────────────────────────────────────────────────────────────────
set -euo pipefail

THREADS="${1:-${THREADS:-8}}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CFG="$ROOT/config/config.yaml"
OUTDIR="$ROOT/results/11_md"
PDB="$ROOT/results/10_structure/construct_refined.pdb"

if ! command -v gmx >/dev/null 2>&1; then
    echo "ERRO: 'gmx' não encontrado no PATH." >&2
    echo "Instale o GROMACS (conda install -c conda-forge gromacs) e rode de novo." >&2
    echo "Nenhum resultado de MD foi gerado — não há como simular sem o motor." >&2
    exit 1
fi

if [[ ! -f "$PDB" ]]; then
    echo "ERRO: modelo do construto ausente: $PDB" >&2
    echo "Rode antes: python scripts/10_structure.py" >&2
    exit 1
fi

# Lê ns e réplicas do config sem depender de parser YAML externo.
NS="$(grep -E '^\s+ns:' "$CFG" | head -1 | tr -dc '0-9')"
REPS="$(grep -E '^\s+replicates:' "$CFG" | head -1 | tr -dc '0-9')"
NS="${NS:-100}"
REPS="${REPS:-3}"
# 2 fs por passo -> nsteps = ns * 500000
NSTEPS=$(( NS * 500000 ))

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "[md] construto=$PDB  ns=$NS  réplicas=$REPS  threads=$THREADS  nsteps=$NSTEPS"

# ── .mdp gerados inline: manter os parâmetros junto do script evita a dessincronia
#    clássica entre o .mdp versionado e o que foi realmente rodado ────────────────
cat > em.mdp <<'EOF'
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.2
rvdw            = 1.2
pbc             = xyz
EOF

cat > nvt.mdp <<'EOF'
integrator      = md
dt              = 0.002
nsteps          = 50000          ; 100 ps
continuation    = no
constraint_algorithm = lincs
constraints     = h-bonds
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.2
rvdw            = 1.2
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1   0.1
ref_t           = 310   310      ; temperatura fisiológica, não 300 K
pcoupl          = no
pbc             = xyz
gen_vel         = yes
gen_temp        = 310
gen_seed        = SEED_PLACEHOLDER
EOF

cat > npt.mdp <<'EOF'
integrator      = md
dt              = 0.002
nsteps          = 50000          ; 100 ps
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.2
rvdw            = 1.2
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1   0.1
ref_t           = 310   310
pcoupl          = C-rescale
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5
pbc             = xyz
gen_vel         = no
EOF

cat > md.mdp <<EOF
integrator      = md
dt              = 0.002
nsteps          = $NSTEPS
nstxout-compressed = 5000        ; frame a cada 10 ps
nstenergy       = 5000
nstlog          = 5000
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.2
rvdw            = 1.2
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1   0.1
ref_t           = 310   310
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5
pbc             = xyz
gen_vel         = no
EOF

# ── Preparo do sistema (uma vez; as réplicas divergem só nas velocidades) ────────
if [[ ! -f solv_ions.gro ]]; then
    # CHARMM36 = opção 1 no menu do pdb2gmx quando o FF está em ./charmm36.ff
    echo "[md] pdb2gmx (CHARMM36 + TIP3P)"
    gmx pdb2gmx -f "$PDB" -o proc.gro -water tip3p -ff charmm36-jul2022 -ignh

    # Dodecaedro: ~29% menos água que um cubo para a mesma distância mínima.
    echo "[md] caixa dodecaédrica, 1.2 nm de folga"
    gmx editconf -f proc.gro -o box.gro -c -d 1.2 -bt dodecahedron

    echo "[md] solvatação"
    gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top

    echo "[md] íons até 0.15 M NaCl (força iônica fisiológica) + neutralização"
    gmx grompp -f em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
    echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top \
        -pname NA -nname CL -neutral -conc 0.15

    echo "[md] minimização de energia"
    gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 2
    gmx mdrun -v -deffnm em -nt "$THREADS"
fi

# ── Réplicas ────────────────────────────────────────────────────────────────────
SUMMARY="$OUTDIR/md_summary.tsv"
printf "replicate\tseed\tns\trmsd_mean_nm\trmsd_sd_nm\trmsf_mean_nm\tgyrate_mean_nm\thbonds_mean\n" > "$SUMMARY"

for ((r=1; r<=REPS; r++)); do
    SEED=$(( 20250719 + r * 7919 ))
    RD="rep${r}"
    mkdir -p "$RD"
    echo "[md] ── réplica $r (seed=$SEED) ──"

    sed "s/SEED_PLACEHOLDER/$SEED/" nvt.mdp > "$RD/nvt.mdp"

    gmx grompp -f "$RD/nvt.mdp" -c em.gro -r em.gro -p topol.top -o "$RD/nvt.tpr" -maxwarn 2
    gmx mdrun -deffnm "$RD/nvt" -nt "$THREADS"

    gmx grompp -f npt.mdp -c "$RD/nvt.gro" -r "$RD/nvt.gro" -t "$RD/nvt.cpt" \
        -p topol.top -o "$RD/npt.tpr" -maxwarn 2
    gmx mdrun -deffnm "$RD/npt" -nt "$THREADS"

    gmx grompp -f md.mdp -c "$RD/npt.gro" -t "$RD/npt.cpt" -p topol.top \
        -o "$RD/md.tpr" -maxwarn 2
    gmx mdrun -deffnm "$RD/md" -nt "$THREADS"

    # Remove saltos de PBC antes de qualquer análise; sem isso o RMSD explode.
    echo -e "Protein\nSystem" | gmx trjconv -s "$RD/md.tpr" -f "$RD/md.xtc" \
        -o "$RD/md_noPBC.xtc" -pbc mol -center

    echo -e "Backbone\nBackbone" | gmx rms -s "$RD/md.tpr" -f "$RD/md_noPBC.xtc" \
        -o "$RD/rmsd.xvg" -tu ns
    echo "Backbone" | gmx rmsf -s "$RD/md.tpr" -f "$RD/md_noPBC.xtc" \
        -o "$RD/rmsf.xvg" -res
    echo "Protein" | gmx gyrate -s "$RD/md.tpr" -f "$RD/md_noPBC.xtc" -o "$RD/gyrate.xvg"
    echo -e "Protein\nProtein" | gmx hbond -s "$RD/md.tpr" -f "$RD/md_noPBC.xtc" \
        -num "$RD/hbond.xvg"

    # Médias das colunas de dados dos .xvg (ignora cabeçalhos @ e #).
    stat_col() {  # $1=arquivo $2=coluna $3=stat(mean|sd)
        awk -v c="$2" -v s="$3" '!/^[@#]/ && NF>=c {
            n++; x=$c; sum+=x; sq+=x*x
        } END {
            if (n==0) { print "NA"; exit }
            m=sum/n;
            if (s=="sd") { v=sq/n-m*m; print (v>0 ? sqrt(v) : 0) } else { print m }
        }' "$1"
    }

    printf "%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\n" \
        "$r" "$SEED" "$NS" \
        "$(stat_col "$RD/rmsd.xvg" 2 mean)" \
        "$(stat_col "$RD/rmsd.xvg" 2 sd)" \
        "$(stat_col "$RD/rmsf.xvg" 2 mean)" \
        "$(stat_col "$RD/gyrate.xvg" 2 mean)" \
        "$(stat_col "$RD/hbond.xvg" 2 mean)" >> "$SUMMARY"
done

echo "[md] concluído. Resumo em $SUMMARY"
cat "$SUMMARY"
