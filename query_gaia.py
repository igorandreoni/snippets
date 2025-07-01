import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table, vstack
from astroquery.gaia import Gaia
import warnings
from tqdm import tqdm
import time
import argparse

def query_gaia_around_coordinates(ra_array, dec_array, radius=2*u.arcsec, 
                                 data_release='DR3', max_sources_per_query=1000,
                                 include_proper_motion=False, mag_limit=None):
    """
    Query Gaia catalog around input coordinates and return matched sources.
    
    Parameters:
    -----------
    ra_array : array-like
        Right Ascension coordinates in degrees
    dec_array : array-like  
        Declination coordinates in degrees
    radius : astropy.units.Quantity, optional
        Search radius around each coordinate (default: 2 arcsec)
    data_release : str, optional
        Gaia data release to query ('DR3', 'EDR3', 'DR2')
    max_sources_per_query : int, optional
        Maximum sources to return per coordinate
    include_proper_motion : bool, optional
        Include proper motion data in results
    mag_limit : float, optional
        Magnitude limit for G-band (e.g., 20.0)
    
    Returns:
    --------
    astropy.table.Table
        Combined results with Gaia sources and input coordinate indices
    """
    
    # Convert inputs to numpy arrays
    ra_array = np.asarray(ra_array)
    dec_array = np.asarray(dec_array)
    
    if len(ra_array) != len(dec_array):
        raise ValueError("RA and Dec arrays must have the same length")
    
    # Set Gaia data release
    if data_release == 'DR3':
        table_name = 'gaiadr3.gaia_source'
    elif data_release == 'EDR3':
        table_name = 'gaiaedr3.gaia_source'
    elif data_release == 'DR2':
        table_name = 'gaiadr2.gaia_source'
    else:
        raise ValueError("data_release must be 'DR3', 'EDR3', or 'DR2'")
    
    print(f"Querying {data_release} for {len(ra_array)} input coordinates")
    print(f"Search radius: {radius}")
    print(f"Maximum sources per query: {max_sources_per_query}")
    
    all_results = []
    failed_queries = []
    
    # Query each coordinate
    for i, (ra, dec) in enumerate(tqdm(zip(ra_array, dec_array), 
                                      total=len(ra_array), 
                                      desc="Querying Gaia")):
        
        try:
            # Create coordinate object
            coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
            
            # Build ADQL query
            radius_deg = radius.to(u.degree).value
            
            # Base columns to select
            columns = [
                'source_id', 'ra', 'dec', 'ra_error', 'dec_error',
                'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag',
                'parallax', 'parallax_error'
            ]
            
            # Add proper motion if requested
            if include_proper_motion:
                columns.extend(['pmra', 'pmdec', 'pmra_error', 'pmdec_error'])
            
            # Build WHERE clause
            where_conditions = [
                f"CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius_deg})) = 1"
            ]
            
            if mag_limit is not None:
                where_conditions.append(f"phot_g_mean_mag < {mag_limit}")
            
            # Construct full query
            query = f"""
            SELECT TOP {max_sources_per_query}
                {', '.join(columns)}
            FROM {table_name}
            WHERE {' AND '.join(where_conditions)}
            ORDER BY phot_g_mean_mag ASC
            """
            
            # Execute query
            job = Gaia.launch_job(query)
            result = job.get_results()
            if len(result) > 0:
                # Add input coordinate index and original coordinates
                result['input_index'] = i
                result['input_ra'] = ra
                result['input_dec'] = dec
                
                # Calculate actual separations
                input_coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
                gaia_coords = SkyCoord(ra=result['ra'], 
                                     dec=result['dec'])
                separations = input_coord.separation(gaia_coords)
                result['separation_arcsec'] = separations.to(u.arcsec).value
                
                all_results.append(result)
            
            # Small delay to be respectful to the server
            time.sleep(0.1)
            
        except Exception as e:
            print(f"Failed to query coordinate {i} (RA={ra:.6f}, Dec={dec:.6f}): {e}")
            failed_queries.append(i)
            continue
    
    # Combine all results
    if all_results:
        combined_results = vstack(all_results)
        print(f"\nQuery completed successfully!")
        print(f"Total Gaia sources found: {len(combined_results)}")
        print(f"Input coordinates with matches: {len(np.unique(combined_results['input_index']))}")
        print(f"Failed queries: {len(failed_queries)}")
        
        return combined_results
    else:
        print("No sources found for any input coordinates")
        return Table()

def analyze_gaia_results(gaia_results):
    """
    Analyze the Gaia query results and provide summary statistics.
    """
    if len(gaia_results) == 0:
        print("No results to analyze")
        return
    
    print("\n=== GAIA QUERY ANALYSIS ===")
    
    # Basic statistics
    print(f"Total sources: {len(gaia_results)}")
    print(f"Unique input coordinates matched: {len(np.unique(gaia_results['input_index']))}")
    
    # Separation statistics
    separations = gaia_results['separation_arcsec']
    print(f"\nSeparation statistics:")
    print(f"  Mean: {np.mean(separations):.3f} arcsec")
    print(f"  Median: {np.median(separations):.3f} arcsec")
    print(f"  Min: {np.min(separations):.3f} arcsec")
    print(f"  Max: {np.max(separations):.3f} arcsec")
    
    # Magnitude statistics
    g_mags = gaia_results['phot_g_mean_mag']
    valid_g = ~np.isnan(g_mags)
    if np.any(valid_g):
        print(f"\nG magnitude statistics:")
        print(f"  Mean: {np.mean(g_mags[valid_g]):.2f}")
        print(f"  Median: {np.median(g_mags[valid_g]):.2f}")
        print(f"  Range: {np.min(g_mags[valid_g]):.2f} - {np.max(g_mags[valid_g]):.2f}")
    
    # Parallax statistics
    parallaxes = gaia_results['parallax']
    valid_plx = ~np.isnan(parallaxes)
    if np.any(valid_plx):
        print(f"\nParallax statistics:")
        print(f"  Sources with parallax: {np.sum(valid_plx)} ({100*np.sum(valid_plx)/len(gaia_results):.1f}%)")
        print(f"  Mean parallax: {np.mean(parallaxes[valid_plx]):.3f} mas")
        print(f"  Median parallax: {np.median(parallaxes[valid_plx]):.3f} mas")
        
        # Distance estimates (for positive parallaxes)
        positive_plx = parallaxes > 0
        if np.any(positive_plx):
            distances = 1000.0 / parallaxes[positive_plx]  # distance in pc
            print(f"  Distance range: {np.min(distances):.1f} - {np.max(distances):.1f} pc")
    
    # Multiple matches per input coordinate
    input_counts = np.bincount(gaia_results['input_index'])
    multiple_matches = np.sum(input_counts > 1)
    if multiple_matches > 0:
        print(f"\nInput coordinates with multiple Gaia matches: {multiple_matches}")
        max_matches = np.max(input_counts)
        print(f"Maximum matches for single coordinate: {max_matches}")

# Example usage
if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Query Gaia catalog around specified coordinates.")
    parser.add_argument('--ra', type=float, nargs='+', required=True, 
                        help='Right Ascension coordinates in degrees (space-separated)')
    parser.add_argument('--dec', type=float, nargs='+', required=True, 
                        help='Declination coordinates in degrees (space-separated)')
    parser.add_argument('--radius', type=float, default=2.0, 
                        help='Search radius in arcseconds (default: 2.0)')
    parser.add_argument('--data_release', type=str, choices=['DR3', 'EDR3', 'DR2'], default='DR3', 
                        help='Gaia data release to query (default: DR3)')
    parser.add_argument('--max_sources', type=int, default=1000, 
                        help='Maximum sources to return per coordinate (default: 1000)')
    parser.add_argument('--include_pm', action='store_true', 
                        help='Include proper motion data in results',
                        default=True)
    parser.add_argument('--mag_limit', type=float, 
                        help='Magnitude limit for G-band (e.g., 20.0)',
                        default=None)

    args = parser.parse_args()
    
    # Query Gaia
    gaia_sources = query_gaia_around_coordinates(
        args.ra, args.dec,
        radius=args.radius * u.arcsec,
        data_release=args.data_release,
        max_sources_per_query=args.max_sources,
        include_proper_motion=args.include_pm,
        mag_limit=args.mag_limit
    )


    print(gaia_sources)
    