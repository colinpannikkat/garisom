"""
Constants for the GARISOM-Torch model.

All physical constants, model parameters, and configuration values.
"""

import torch
from typing import Final

# Device configuration - prefer MPS (Apple Silicon) > CUDA > CPU
def _get_best_device() -> torch.device:
    """Determine the best available device."""
    if torch.backends.mps.is_available():
        return torch.device("mps")
    elif torch.cuda.is_available():
        return torch.device("cuda")
    else:
        return torch.device("cpu")

DEVICE: torch.device = _get_best_device()
# Note: MPS has limited float64 support, use float32 for MPS
DTYPE: torch.dtype = torch.float32 if DEVICE.type == "mps" else torch.float64

print(f"GARISOM-Torch using device: {DEVICE} with dtype: {DTYPE}")

# Mathematical constants
PI: Final[float] = 3.14159265358979323846

# Physical constants
SBC: Final[float] = 5.67e-8          # Stefan-Boltzmann constant (W m^-2 K^-4)
SHA: Final[float] = 29.3             # Specific heat of air (J mol^-1 C^-1)
GAS: Final[float] = 8.3144598        # Universal gas constant (J mol^-1 K^-1)
OA: Final[float] = 0.21              # Mole fraction of O2
SOLAR: Final[float] = 1362.0         # Solar constant (W m^-2)

# Leaf optical properties
ABS_SOLAR: Final[float] = 0.5        # Absorptivity of solar for leaves
ABS_PAR: Final[float] = 0.8          # Absorptivity of PAR for leaves
ABS_NIR: Final[float] = 0.2          # Absorptivity of near infrared for leaves

# Model thresholds and limits
GMAX: Final[float] = 1e6             # Maximum G, wet soil, vpd=0 (kg m^-2 hr^-1 basal area)
MIN_WIND_THRESH: Final[float] = 0.4515  # Minimum wind threshold (m s^-1)
DPA_MAX_CUTOFF: Final[float] = 1.1   # Cutoff for stopping dpamax search
PROFT_MAX_RUN_MEAN: Final[int] = 1   # Running mean for profit maximization

# Numerical integration parameters
TRAP_ITER_MAX: Final[int] = 70       # Max iterations for trapezoid rule
FLOW_ITER_LIMIT: Final[int] = 100000 # Max iterations for flow rate calculation
XYLEM_EPSX: Final[float] = 0.0001    # Acceptable error for xylem E integral
RHIZO_EPSX: Final[float] = 0.001     # Acceptable error for rhizosphere E integral

# Array sizes
CURVE_MAX: Final[int] = 10001        # Reduced from 100001 for memory efficiency
MAX_YEARS: Final[int] = 90
MAX_LAYERS: Final[int] = 10          # Maximum number of soil layers

# Stage IDs for model configuration
STAGE_ID_NONE: Final[int] = 0
STAGE_ID_HIST_OPT: Final[int] = 1
STAGE_ID_HIST_STRESS: Final[int] = 2
STAGE_ID_FUT_OPT: Final[int] = 3
STAGE_ID_FUT_STRESS: Final[int] = 4
STAGE_ID_FUT_STRESS_NOACCLIM: Final[int] = 5

# Default atmospheric parameters
DEFAULT_PATM: Final[float] = 101.3   # Atmospheric pressure (kPa)

# Van Genuchten soil parameters by texture type
# Format: {texture: (alpha [MPa^-1], n, ksat [kg hr^-1 MPa^-1 m^-2], theta_sat)}
VAN_GENUCHTEN_PARAMS: dict = {
    "sand": (1479.5945, 2.68, 30305.88, 0.43),
    "loamy sand": (1265.3084, 2.28, 14897.84, 0.41),
    "sandy loam": (765.3075, 1.89, 4510.168, 0.41),
    "loam": (367.3476, 1.56, 1061.216, 0.43),
    "silt": (163.2656, 1.37, 255.1, 0.46),
    "silt loam": (204.082, 1.41, 459.18, 0.45),
    "sandy clay loam": (602.0419, 1.48, 1336.724, 0.39),
    "clay loam": (193.8779, 1.31, 265.304, 0.41),
    "silty clay loam": (102.041, 1.23, 71.428, 0.43),
    "sandy clay": (275.5107, 1.23, 122.448, 0.38),
    "silty clay": (51.0205, 1.09, 20.408, 0.36),
    "clay": (81.6328, 1.09, 204.08, 0.38),
}


def get_device() -> torch.device:
    """Get the current device for tensor allocation."""
    return DEVICE


def get_dtype() -> torch.dtype:
    """Get the current dtype for tensor creation."""
    return DTYPE


def tensor(data, requires_grad: bool = False) -> torch.Tensor:
    """Create a tensor with consistent device and dtype settings.
    
    Args:
        data: Input data (scalar, list, numpy array, or tensor)
        requires_grad: Whether to enable gradient computation
        
    Returns:
        Tensor on the configured device with configured dtype
    """
    if isinstance(data, torch.Tensor):
        return data.to(device=DEVICE, dtype=DTYPE).requires_grad_(requires_grad)
    return torch.tensor(data, device=DEVICE, dtype=DTYPE, requires_grad=requires_grad)


def zeros(*size, requires_grad: bool = False) -> torch.Tensor:
    """Create a zero tensor with consistent device and dtype."""
    return torch.zeros(*size, device=DEVICE, dtype=DTYPE, requires_grad=requires_grad)


def ones(*size, requires_grad: bool = False) -> torch.Tensor:
    """Create a ones tensor with consistent device and dtype."""
    return torch.ones(*size, device=DEVICE, dtype=DTYPE, requires_grad=requires_grad)


def empty(*size, requires_grad: bool = False) -> torch.Tensor:
    """Create an empty tensor with consistent device and dtype."""
    return torch.empty(*size, device=DEVICE, dtype=DTYPE, requires_grad=requires_grad)
