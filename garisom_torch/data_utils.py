"""
Data Loading Utilities for GARISOM-Torch.

This module provides functions for loading:
    - Model parameters from CSV
    - Climate forcing data from CSV
    - Configuration settings from CSV
    
The utilities map the original C++ column naming conventions
to the PyTorch model parameter structure.
"""

import torch
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Optional, Union, Tuple
from dataclasses import dataclass

from .constants import PI, tensor, get_device, get_dtype


# Parameter name mapping from CSV columns to model attributes
PARAM_MAPPING = {
    # Site parameters
    'i_latitude': 'lat',
    'i_longitude': 'longitude',
    'i_slopeI': 'slope',
    'i_slopeA': 'slope_asp',
    'i_elevation': 'alt',
    
    # Atmosphere
    'i_atmTrans': 'tau',
    'i_solarNoon': 'tsn_corr',
    'i_co2AmbPPM': 'ca',
    'i_emiss': 'emiss',
    
    # Soil
    'i_layers': 'layers',
    'i_rockFrac': 'rock_frac',
    'i_rhizoPer': 'rhizo_targ',
    'i_fieldCapPercInit': 'ffc',
    'i_fieldCapFrac': 'field_cap_frac',
    'i_gWaterP': 'p_ground',
    'i_gWaterDist': 'ground_distance',
    'i_soilAbsSol': 'soil_abs_sol',
    
    # Stand parameters
    'i_treeToPhotoLAI': 'tree_to_photo_lai',
    'i_leafAreaIndex': 'lai',
    'i_leafAngleParam': 'leaf_angle_param',
    'i_leafPerBasal': 'leaf_per_basal',
    'i_height': 'height',
    'i_aspect': 'aspect',
    'i_rootDepth': 'root_depth',
    'i_rootBeta': 'beta',
    'i_leafWidth': 'leaf_width',
    'i_baperga': 'ba_per_ga',
    'i_soilXHeight': 'xh',
    
    # Hydraulics
    'i_kmaxTree': 'ksatp',
    'i_pinc': 'p_inc',
    'i_leafPercRes': 'leaf_percent_res',
    
    # Weibull parameters
    'i_cr': 'c_root',
    'i_br': 'b_root',
    'i_cs': 'c_stem',
    'i_bs': 'b_stem',
    'i_cl': 'c_leaf',
    'i_bl': 'b_leaf',
    
    # Photosynthesis
    'i_lightComp': 'light_comp',
    'i_qMax': 'q_max',
    'i_vmax25': 'v_max25',
    'i_jmax25': 'j_max25',
    'i_kc25': 'kc_25',
    'i_ko25': 'ko_25',
    'i_comp25': 'comp_25',
    'i_thetaC': 'theta_c',
    'i_havmax': 'ha_vmax',
    'i_hdvmax': 'hd_vmax',
    'i_svvmax': 'sv_vmax',
    'i_hajmax': 'ha_jmax',
    'i_hdjmax': 'hd_jmax',
    'i_svjmax': 'sv_jmax',
    'i_lightCurv': 'light_curv',
    
    # Additional hydraulic parameters
    'i_sapwoodT': 'sapwood_t',
    'i_conduitDiam': 'conduit_d',
}

# Climate data column mapping
CLIMATE_MAPPING = {
    'Year': 'year',
    'Day': 'julian_day',
    'Hour': 'hour',
    'Solar_Wm2': 'solar',
    'Rain_mm': 'rain',
    'Wind_ms.1': 'wind',
    'Tair_C': 'tair',
    'Tsoil_C': 'tsoil',
    'D_kPa': 'vpd',
}


@dataclass
class ClimateData:
    """
    Container for climate forcing data as PyTorch tensors.
    
    Attributes:
        year: Year values
        julian_day: Julian day of year (1-365)
        hour: Hour of day (0-23)
        solar: Solar radiation (W/m²)
        rain: Rainfall (mm)
        wind: Wind speed (m/s)
        tair: Air temperature (°C)
        tsoil: Soil temperature (°C)
        vpd: Vapor pressure deficit (kPa)
        n_timesteps: Number of timesteps
    """
    year: torch.Tensor
    julian_day: torch.Tensor
    hour: torch.Tensor
    solar: torch.Tensor
    rain: torch.Tensor
    wind: torch.Tensor
    tair: torch.Tensor
    tsoil: torch.Tensor
    vpd: torch.Tensor
    
    @property
    def n_timesteps(self) -> int:
        return len(self.year)
    
    def __len__(self) -> int:
        return self.n_timesteps
    
    def __getitem__(self, idx: int) -> Dict[str, torch.Tensor]:
        """Get single timestep as dictionary."""
        return {
            'year': self.year[idx],
            'julian_day': self.julian_day[idx],
            'hour': self.hour[idx],
            'solar': self.solar[idx],
            'rain': self.rain[idx],
            'wind': self.wind[idx],
            'tair': self.tair[idx],
            'tsoil': self.tsoil[idx],
            'vpd': self.vpd[idx],
        }
    
    def to(self, device: torch.device) -> 'ClimateData':
        """Move all tensors to device."""
        return ClimateData(
            year=self.year.to(device),
            julian_day=self.julian_day.to(device),
            hour=self.hour.to(device),
            solar=self.solar.to(device),
            rain=self.rain.to(device),
            wind=self.wind.to(device),
            tair=self.tair.to(device),
            tsoil=self.tsoil.to(device),
            vpd=self.vpd.to(device),
        )


@dataclass
class ModelConfig:
    """
    Container for model configuration settings.
    
    Attributes:
        ground_water_enabled: Enable groundwater flow
        soil_redistribution_enabled: Enable soil water redistribution
        soil_evaporation_enabled: Enable soil evaporation
        rain_enabled: Enable rainfall inputs
        refilling_enabled: Enable xylem refilling
        predawns_mode: Use measured predawn potentials
        cavitation_fatigue: Enable cavitation fatigue/hysteresis
        stem_only: Fatigue in stems only (not roots)
        multiple_species: Enable multi-species mode
        species_number: Number of species
        climate_data_path: Path to climate forcing data
        growing_season_path: Path to growing season limits
    """
    ground_water_enabled: bool = False
    soil_redistribution_enabled: bool = True
    soil_evaporation_enabled: bool = True
    rain_enabled: bool = True
    refilling_enabled: bool = False
    predawns_mode: bool = False
    cavitation_fatigue: bool = False
    stem_only: bool = True
    multiple_species: bool = False
    species_number: int = 1
    climate_data_path: str = ""
    growing_season_path: str = ""


def load_parameters(filepath: Union[str, Path], 
                    species_index: int = 0) -> Dict[str, float]:
    """
    Load model parameters from CSV file.
    
    Args:
        filepath: Path to parameters CSV file
        species_index: Row index for species (0-based)
        
    Returns:
        Dictionary of parameter names and values
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"Parameter file not found: {filepath}")
    
    # Read CSV
    df = pd.read_csv(filepath)
    
    # Get row for species
    if species_index >= len(df):
        raise IndexError(f"Species index {species_index} out of range (max {len(df)-1})")
    
    row = df.iloc[species_index]
    
    # Map to model parameters
    params = {}
    
    for csv_col, model_param in PARAM_MAPPING.items():
        if csv_col in row.index:
            value = row[csv_col]
            
            # Handle special transformations
            if csv_col == 'i_latitude':
                value = PI / 180 * value  # Convert to radians
            elif csv_col == 'i_co2AmbPPM':
                pass  # Keep as ppm, model will convert
            elif csv_col == 'i_rhizoPer':
                value = value / 100.0  # Convert percent to fraction
            elif csv_col == 'i_fieldCapPercInit':
                value = value / 100.0  # Convert percent to fraction
            elif csv_col == 'i_baperga':
                value = value * 0.0001  # Convert m²/ha to m²/m²
            elif csv_col == 'i_leafWidth':
                value = value * 0.72  # Apply characteristic dimension factor
            elif csv_col == 'i_soilXHeight':
                value = value * 100.0  # Convert to cm
                
            params[model_param] = float(value)
    
    # Add texture if present
    if 'i_texture' in row.index:
        params['texture'] = str(row['i_texture'])
    
    # Add species info if present
    if 'i_sp' in row.index:
        params['species'] = str(row['i_sp'])
    if 'i_site' in row.index:
        params['site'] = str(row['i_site'])
        
    return params


def load_climate_data(filepath: Union[str, Path],
                      start_row: Optional[int] = None,
                      end_row: Optional[int] = None) -> ClimateData:
    """
    Load climate forcing data from CSV file.
    
    Args:
        filepath: Path to climate data CSV file
        start_row: Starting row index (optional)
        end_row: Ending row index (optional)
        
    Returns:
        ClimateData object with forcing tensors
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"Climate file not found: {filepath}")
    
    # Read CSV
    df = pd.read_csv(filepath)
    
    # Subset if requested
    if start_row is not None or end_row is not None:
        start = start_row or 0
        end = end_row or len(df)
        df = df.iloc[start:end]
    
    # Map columns (handle various naming conventions)
    def get_column(df, possible_names):
        for name in possible_names:
            if name in df.columns:
                return df[name].values
        raise KeyError(f"Could not find column matching: {possible_names}")
    
    # Extract data with flexible column naming
    year = get_column(df, ['Year', 'year', 'YEAR'])
    julian_day = get_column(df, ['Day', 'julian-day', 'Julian_Day', 'DOY', 'doy'])
    hour = get_column(df, ['Hour', 'hour', 'standard-time', 'Time'])
    solar = get_column(df, ['Solar_Wm2', 'solar', 'Solar', 'SOLAR', 'Radiation'])
    rain = get_column(df, ['Rain_mm', 'rain', 'Rain', 'RAIN', 'Precip'])
    wind = get_column(df, ['Wind_ms.1', 'wind', 'Wind', 'WIND', 'Wind_ms'])
    tair = get_column(df, ['Tair_C', 'T-air', 'Tair', 'TAIR', 'AirTemp'])
    tsoil = get_column(df, ['Tsoil_C', 'T-soil', 'Tsoil', 'TSOIL', 'SoilTemp'])
    vpd = get_column(df, ['D_kPa', 'D-MD', 'VPD', 'vpd', 'D'])
    
    # Convert to tensors
    device = get_device()
    dtype = get_dtype()
    
    return ClimateData(
        year=torch.tensor(year, device=device, dtype=torch.long),
        julian_day=torch.tensor(julian_day, device=device, dtype=torch.long),
        hour=torch.tensor(hour, device=device, dtype=dtype),
        solar=torch.tensor(solar, device=device, dtype=dtype),
        rain=torch.tensor(rain, device=device, dtype=dtype),
        wind=torch.tensor(wind, device=device, dtype=dtype),
        tair=torch.tensor(tair, device=device, dtype=dtype),
        tsoil=torch.tensor(tsoil, device=device, dtype=dtype),
        vpd=torch.tensor(vpd, device=device, dtype=dtype),
    )


def load_configuration(filepath: Union[str, Path],
                       config_index: int = 0) -> ModelConfig:
    """
    Load model configuration from CSV file.
    
    Args:
        filepath: Path to configuration CSV file
        config_index: Row index for configuration (0-based)
        
    Returns:
        ModelConfig object with settings
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"Configuration file not found: {filepath}")
    
    # Read CSV
    df = pd.read_csv(filepath)
    
    if config_index >= len(df):
        raise IndexError(f"Config index {config_index} out of range")
    
    row = df.iloc[config_index]
    
    def parse_bool(val):
        if isinstance(val, bool):
            return val
        return str(val).lower() in ('y', 'yes', 'true', '1', 'on')
    
    config = ModelConfig(
        ground_water_enabled=parse_bool(row.get('i_gWaterEnable', 'n')),
        soil_redistribution_enabled=parse_bool(row.get('i_soilRedEnable', 'y')),
        soil_evaporation_enabled=parse_bool(row.get('i_soilEvapEnable', 'y')),
        rain_enabled=parse_bool(row.get('i_rainEnable', 'y')),
        refilling_enabled=parse_bool(row.get('i_refilling', 'n')),
        predawns_mode=parse_bool(row.get('i_predawnsMode', 'n')),
        cavitation_fatigue=parse_bool(row.get('i_cavitFatigue', 'n')),
        stem_only=parse_bool(row.get('i_stemOnly', 'y')),
        multiple_species=parse_bool(row.get('i_multipleSP', 'n')),
        species_number=int(row.get('i_speciesN', 1)),
        climate_data_path=str(row.get('i_ClimateData', '')),
        growing_season_path=str(row.get('i_GSData', '')),
    )
    
    return config


def create_data_loader(climate_data: ClimateData,
                       batch_size: int = 1,
                       shuffle: bool = False) -> torch.utils.data.DataLoader:
    """
    Create a PyTorch DataLoader for climate data.
    
    Args:
        climate_data: ClimateData object
        batch_size: Batch size for loading
        shuffle: Whether to shuffle data
        
    Returns:
        DataLoader for iterating over timesteps
    """
    # Stack all tensors
    data = torch.stack([
        climate_data.year.float(),
        climate_data.julian_day.float(),
        climate_data.hour,
        climate_data.solar,
        climate_data.rain,
        climate_data.wind,
        climate_data.tair,
        climate_data.tsoil,
        climate_data.vpd,
    ], dim=1)
    
    dataset = torch.utils.data.TensorDataset(data)
    
    return torch.utils.data.DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
    )


def save_results(results: Dict[str, torch.Tensor],
                 filepath: Union[str, Path],
                 climate_data: Optional[ClimateData] = None):
    """
    Save model results to CSV file.
    
    Args:
        results: Dictionary of result tensors (timestep, value)
        filepath: Output file path
        climate_data: Optional climate data to include timestamps
    """
    filepath = Path(filepath)
    
    # Convert tensors to numpy
    data = {}
    for key, val in results.items():
        if isinstance(val, torch.Tensor):
            data[key] = val.detach().cpu().numpy()
        else:
            data[key] = val
            
    # Add climate timestamps if available
    if climate_data is not None:
        data['year'] = climate_data.year.cpu().numpy()
        data['julian_day'] = climate_data.julian_day.cpu().numpy()
        data['hour'] = climate_data.hour.cpu().numpy()
        
    # Create DataFrame and save
    df = pd.DataFrame(data)
    df.to_csv(filepath, index=False)
    print(f"Results saved to {filepath}")
