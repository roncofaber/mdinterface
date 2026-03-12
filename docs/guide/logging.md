# Logging

`mdinterface` uses Python's standard `logging` module. By default all log output is suppressed (a `NullHandler` is installed at import time, following standard library practice).

## Enabling output

**Package-wide** — the recommended approach:

```python
import mdinterface

mdinterface.set_verbosity(1)        # INFO  (normal detail)
mdinterface.set_verbosity(2)        # DEBUG (maximum detail)
mdinterface.set_verbosity(0)        # WARNING (quiet)
mdinterface.set_verbosity("DEBUG")  # string form
mdinterface.set_verbosity(True)     # same as 1 / INFO
mdinterface.set_verbosity(False)    # same as 0 / WARNING
```

**Via SimCell constructor** — convenient for one-off scripts:

```python
from mdinterface import SimCell

simbox = SimCell(xysize=[15, 15], verbose=True)   # INFO
simbox = SimCell(xysize=[15, 15], verbose=2)      # DEBUG
```

## Verbosity levels

| Value | Level | Typical output |
|-------|-------|----------------|
| `0` / `False` | WARNING | Only warnings and errors |
| `1` / `True` | INFO | Build summary, layer sizes, molecule counts |
| `2` | DEBUG | PACKMOL details, internal operations |
| `"DEBUG"` / `"INFO"` / ... | raw Python level | Passed directly to `logging` |

## Log format

All messages are prefixed with `[mdi]` and a compact 4-character level name:

```
[mdi] INFO | ===  Build  ================================
[mdi] INFO |
[mdi] INFO | --  Layer [1/3]  --------------------------
[mdi] INFO |   >> Au (111):  14.421 x 14.421 x 6.657 Å,  480 atoms
[mdi] INFO |   >> layer z: +6.66 Å  |  total z: 6.66 Å
```

## Integration with existing logging config

Since `mdinterface` uses a dedicated `"mdinterface"` logger subtree, it coexists cleanly with any logging configuration you already have in your application. `set_verbosity` only affects the `"mdinterface"` logger family and does not touch the root logger.
