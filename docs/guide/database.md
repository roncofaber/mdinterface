# Database

`mdinterface` ships a built-in database of common species with pre-defined force-field parameters. All entries are importable from `mdinterface.database`.

## Metals

```python
from mdinterface.database import Metal111

gold    = Metal111("Au")
silver  = Metal111("Ag")
platinum = Metal111("Pt")
copper  = Metal111("Cu")
```

`Metal111` generates an FCC (111) surface slab. The element symbol selects the lattice parameter and Lennard-Jones parameters.

## Water models

```python
from mdinterface.database import Water

water_spce  = Water(model="ewald") # modified tip3p model
water_tip4p = Water(model="spce")  # SPC/E water
```

## Ions

```python
from mdinterface.database import Ion

na = Ion("Na", ffield="Cheatham")
cl = Ion("Cl", ffield="Cheatham")
li = Ion("Li", ffield="Cheatham")
k  = Ion("K",  ffield="Cheatham")
```

Currently, the following force field parameters for monovalent ions have been implemented:

- **Aqvist**   : J. Phys. Chem. B 2008, [https://pubs.acs.org/doi/10.1021/jp8001614](https://pubs.acs.org/doi/10.1021/jp8001614),
- **Jorgensen**: J. Chem. Theory Comput. 2006, [https://pubs.acs.org/doi/10.1021/ct600252r](https://pubs.acs.org/doi/10.1021/ct600252r)
- **Cheatham** : J. Phys. Chem. B 2008, [https://pubs.acs.org/doi/10.1021/jp8001614](https://pubs.acs.org/doi/10.1021/jp8001614)
- **Sengupta** : J. Chem. Inf. Model. 2021, [https://pubs.acs.org/doi/10.1021/acs.jcim.0c01390](https://pubs.acs.org/doi/10.1021/acs.jcim.0c01390)
- **Dang**     : J. Chem. Phys. 1992/1994
- **OPLS-AA**  : J. Chem. Theory Comput. 2009, [https://pubs.acs.org/doi/10.1021/ct900009a](https://pubs.acs.org/doi/10.1021/ct900009a)

Special ion species also available:

```python
from mdinterface.database import Perchlorate, Hydronium, Hydroxide
```

If the parameter set you are looking for are not present, you can always create a Specie explicitly:

```python
from mdinterface import Specie

Li = Specie("Li", charges=1, lj={"Li": [0.33673, 1.40940]})
```



## Noble gases

```python
from mdinterface.database import Neon, Argon, Krypton, Xenon, NobleGas

ar = Argon()
xe = Xenon()
# or generically:
gas = NobleGas("Kr")
```

## Graphene

```python
from mdinterface.database import Graphene

grap = Graphene()
```

## Small molecules

```python
from mdinterface.database import Oxygen, Hydrogen, Nitrogen
```
