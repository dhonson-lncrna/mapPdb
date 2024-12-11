from chimerax.core.commands import run

with open(chain_file,'r') as chf:
    chain_ls = [i.strip() for i in chf.readlines()]

with open(cmap_file,'r') as cmf:
    cmap_ls = [i.strip() for i in cmf.readlines()]

chain_ls = [i.split(',') for i in chain_ls]
cmap_ls = [i.split(',') for i in cmap_ls]

# Open file and select all chains
run(session, f'open '+pdbid)
run(session, f'select all')

# Remove chains for analysis
for ch in chain_ls:
    run(session, f'select subtract {ch[0]}')

# Delete unwanted chains
run(session, f'delete sel')

# Add silhouette
run(session, f'graphics silhouette true width 1')

# Style RNA and protein chains
for ch in chain_ls:
    run(session, f'color {ch[0]} {ch[1]}')
    if ch[0] in rna_chains:
        pass
    else:
        run(session, f'transparency {ch[0]} 30 all')
        run(session, f'hide {ch[0]} atoms')
        run(session, f'show {ch[0]} surfaces')

for r in rna_chains:
    run(session, f'hide {r} ribbons')
    run(session, f'style {r} stick')
    run(session, f'hide {r} atoms')
    run(session, f'hide {r} surfaces')
    run(session, f'renumber {r}')
    #run(session, f'transparency {r} 100 all')

# Color RNA based on SPIDR data
for cm in cmap_ls:
    run(session, f'color {cm[0]} {cm[1]}')
    run(session, f'show {cm[0]} atoms')

# Save session
run(session, f'save {cxs_file}')


