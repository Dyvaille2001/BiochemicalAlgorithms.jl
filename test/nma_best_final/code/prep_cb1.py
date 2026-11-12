# prep_cb1.py   usage: python prep_cb1.py in.pdb out.pdb seed
import sys
import time, os
from openmm import app, unit, LangevinMiddleIntegrator, CustomExternalForce
from pdbfixer import PDBFixer

inp, out, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])

# 1) keep only the CB1 residues; the 306/332 gap becomes chain A / chain B
keep = []
for l in open(inp):
    if l.startswith("ATOM"):
        r = int(l[22:26])
        if 99 <= r <= 306:   keep.append(l[:21] + "A" + l[22:])
        elif 332 <= r <= 414: keep.append(l[:21] + "B" + l[22:])
open("_tmp_in.pdb", "w").writelines(keep + ["END\n"])

# 2) rebuild missing atoms, add H (no loop building across the gap)
fx = PDBFixer("_tmp_in.pdb")
fx.findMissingResidues(); fx.missingResidues = {}
fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
fx.findMissingAtoms(); fx.addMissingAtoms()
fx.addMissingHydrogens(7.4)

# 3) force field with backbone restraints
ff = app.ForceField("amber14-all.xml", "implicit/obc2.xml")
system = ff.createSystem(fx.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
rest = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
rest.addGlobalParameter("k", 1000.0)            # kJ/mol/nm^2
for n in ("x0", "y0", "z0"): rest.addPerParticleParameter(n)
for a in fx.topology.atoms():
    if a.name in ("N", "CA", "C", "O"):
        p = fx.positions[a.index].value_in_unit(unit.nanometer)
        rest.addParticle(a.index, [p[0], p[1], p[2]])
system.addForce(rest)

# 4) minimise, short MD with this seed, minimise again
#integ = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
#integ.setRandomNumberSeed(seed)
#sim = app.Simulation(fx.topology, system, integ)
#sim.context.setPositions(fx.positions)
#sim.minimizeEnergy(maxIterations=2000)
#sim.context.setVelocitiesToTemperature(300*unit.kelvin, seed)
#sim.step(50000)                                   # 100 ps
#sim.minimizeEnergy(tolerance=1*unit.kilojoule_per_mole/unit.nanometer, maxIterations=5000)

t0 = time.time()
print("minimising 1...", flush=True)
sim.minimizeEnergy(maxIterations=2000)
print(f"  done {time.time()-t0:.0f}s", flush=True)

sim.context.setVelocitiesToTemperature(300*unit.kelvin, seed)
total = int(os.environ.get("MD_STEPS", 50000))
chunk = max(total // 10, 500)
done = 0
while done < total:
    n = min(chunk, total - done)
    sim.step(n); done += n
    print(f"MD {done}/{total} steps, {time.time()-t0:.0f}s", flush=True)

print("minimising 2...", flush=True)
sim.minimizeEnergy(tolerance=1*unit.kilojoule_per_mole/unit.nanometer, maxIterations=5000)
print(f"  done {time.time()-t0:.0f}s", flush=True)

# 5) write heavy atoms only, chain A, original residue numbers
mod = app.Modeller(fx.topology, sim.context.getState(getPositions=True).getPositions())
mod.delete([a for a in mod.topology.atoms() if a.element is not None and a.element.symbol == "H"])
app.PDBFile.writeFile(mod.topology, mod.positions, open("_tmp_out.pdb", "w"), keepIds=True)
with open(out, "w") as o:
    for l in open("_tmp_out.pdb"):
        if l.startswith(("TER", "CONECT", "ANISOU")): continue
        if l.startswith("ATOM"): l = l[:21] + "A" + l[22:]
        o.write(l)
