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

water_spce  = Water(model="ewald")    # SPC/E
water_tip4p = Water(model="tip4p")    # TIP4P/2005
```

## Ions

```python
from mdinterface.database import Ion

na = Ion("Na", ffield="Cheatham")
cl = Ion("Cl", ffield="Cheatham")
li = Ion("Li", ffield="Cheatham")
k  = Ion("K",  ffield="Cheatham")
```

`lookup_parameters` can be used to check what force fields are available for a given ion:

```python
from mdinterface.database import lookup_parameters
lookup_parameters("Na")
```

Special ion species also available:

```python
from mdinterface.database import Perchlorate, Hydronium, Hydroxide
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
