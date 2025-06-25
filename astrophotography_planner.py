import pandas as pd
import requests
import datetime
from timezonefinder import TimezoneFinder
import yaml
import warnings
import pytz
import re

from astroplan import Observer, FixedTarget
from astroplan.constraints import AltitudeConstraint, AtNightConstraint, MoonSeparationConstraint, MoonIlluminationConstraint
from astroplan.utils import time_grid_from_range

from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import NonRotationTransformationWarning
import astropy.units as u

from astroquery.simbad import Simbad

from geopy.geocoders import Nominatim

from ics import Calendar, Event
import matplotlib.pyplot as plt
import calendar
from collections import defaultdict

# supress warnings about slight inaccuracy in Angular separations. Only matters for planets and other non-DSO targets
warnings.filterwarnings('ignore', category=NonRotationTransformationWarning)

# Add necessary fields to Simbad query
Simbad.add_votable_fields('ra(d)', 'dec(d)', 'ids')

def get_config():
    # Open yaml file and read it in to a sict
    with open('config.yaml', 'r') as file:
        config = yaml.safe_load(file)

    # Printing
    print("Program Configuration")
    for key, value in config.items():
        print(f"{key}: {value}")
    print("\n")

    return config

def get_observer(location_name):
    # Get lat and lon from name
    loc = Nominatim(user_agent="GetLoc")
    geo_location = loc.geocode(location_name)
    latitude = geo_location.latitude
    longitude = geo_location.longitude

    # Loop up timezone from lat adn lon
    tz_finder = TimezoneFinder()
    timezone = tz_finder.timezone_at(lat=latitude, lng=longitude)

    # Printing
    print('Observer Location')
    print(f"Latitude: {latitude}")
    print(f"Longitude: {longitude}")
    print(f"Timezone: {timezone}")
    print("\n")

    # Build observer class and return
    observer = Observer(latitude=latitude, longitude=longitude, name=location_name, timezone=timezone)
    
    return observer

def get_targets(config):
    targets = []

    print("Target List:")

    for name in config['targets']:
        result = Simbad.query_object(name)
        if result is None:
            print(f"Could not resolve target: {name}")
            continue
        
        main_id = result['main_id'][0].decode('utf-8') if isinstance(result['main_id'][0], bytes) else result['main_id'][0]
        main_id = re.sub(r'\s+', ' ', main_id).strip().replace('NAME ', '') # Clean multiple white spaces        
        alias = get_pretty_name_from_ids(result['ids'][0], main_id)
        
        if main_id != alias:
            name = f"{alias} [{main_id}]"
        else:
            name = main_id

        ra = result['ra'][0] * u.deg
        dec = result['dec'][0] * u.deg
        coord = SkyCoord(ra=ra, dec=dec)

        # Create FixedTarget with main_id as name
        target = FixedTarget(name=name, coord=coord)
        targets.append(target)

        # Format RA
        ra_str = target.coord.ra.to_string(unit='hourangle', sep=':')

        # Format Dec with degree, arcmin, arcsec
        deg, arcmin, arcsec = target.coord.dec.dms
        sign = '+' if deg >= 0 else '-'
        dec_str = f"{sign}{abs(deg):.0f}° {abs(arcmin):.0f}' {abs(arcsec):.2f}\""

        print(f"Name: {name}")
        print(f"RA: {ra_str}")
        print(f"Dec: {dec_str}\n")

    return targets

def get_pretty_name_from_ids(ids_string, main_id):
    aliases = [alias.strip() for alias in ids_string.split('|')]

    name_aliases = [
        alias.replace('NAME ', '') 
        for alias in aliases 
        if alias.startswith('NAME ') and main_id.replace(' ', '') not in alias.replace(' ', '')
    ]
    if name_aliases:
        return name_aliases[0]
    
    return main_id  # No pretty name found

def make_ics(results):
    df = pd.DataFrame(results)
    df['date'] = pd.to_datetime(df['date'])

    calendar = Calendar()

    for _, row in df.iterrows():
        event = Event()
        event.name = f"Observe {row['target']}"
        event.begin = row['observation_start']
        event.duration = {"hours": int(row['observable_hours'])}
        event.description = f"Target: {row['target']}\nObservation window begins: {row['observation_start']}\nEstimated observable time: {row['observable_hours']:.2f} hours"
        calendar.events.add(event)

    with open("astro_planner_targets.ics", "w") as f:
        f.writelines(calendar)

def make_calendar_images(results, year):
    # Organize results by date
    observations_by_day = defaultdict(list)
    for row in results:
        date_obj = datetime.datetime.strptime(row['date'], "%Y-%m-%d").date()
        observations_by_day[date_obj].append(row)

    # Loop over each month
    for month in range(1, 13):
        cal = calendar.Calendar(firstweekday=0)  # Monday start
        month_days = cal.monthdatescalendar(year, month)

        fig, ax = plt.subplots(figsize=(24, 14))
        ax.set_axis_off()
        ax.set_title(calendar.month_name[month] + f" {year}", fontsize=16, fontweight='bold')

        table_data = []
        for week in month_days:
            week_row = []
            for day in week:
                cell = f"{day.day}"
                if day in observations_by_day:
                    obs_info = "\n".join(
                        f"{obs['target']} - {obs['observation_start'].split()[1]}"
                        for obs in observations_by_day[day]
                    )
                    cell += f"\n{obs_info}"
                week_row.append(cell)
            table_data.append(week_row)

        # Create table
        table = ax.table(cellText=table_data,
                         colLabels=["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"],
                         cellLoc='left',
                         loc='center')
        table.scale(1, 7)

        #Set font size
        for (row, col), cell in table.get_celld().items():
            if row > 0:  # Skip header
                cell.PAD = 0.1
                cell.get_text().set_fontsize(42)
            else:
                cell.PAD = 0.1
                cell.get_text().set_fontsize(42)

        # Bold observation days
        for (i, week) in enumerate(month_days):
            for (j, day) in enumerate(week):
                if day in observations_by_day:
                    table[(i + 1, j)].set_facecolor('#cce5ff')

        plt.tight_layout()
        plt.savefig(f"calendar_{year}_{month:02d}.png", dpi=150)
        plt.close()

def main():
    # Load configuration
    config = get_config()

    # Build observer
    observer = get_observer(config['location_name'])

    # Load target names and resolve coordinates
    targets = get_targets(config)

    # Constraints
    constraints = [
        AltitudeConstraint(min=config['min_altitude_deg'] * u.deg, max=config['max_altitude_deg'] * u.deg),
        AtNightConstraint.twilight_astronomical(),
        MoonSeparationConstraint(min=config['moon_sep_deg'] * u.deg),
        MoonIlluminationConstraint(max=config['max_moon_illum'])
    ]

    # Observation date
    start_date = datetime.date(config['calendar_year'], 1, 1)
    end_date = datetime.date(config['calendar_year'], 12, 31)
    delta = (end_date - start_date).days + 1
    dates = [start_date + datetime.timedelta(days=i) for i in range(delta)]

    # Minimum observable time
    min_minutes = int(config['min_minutes_observable']) * u.minute

    # Results
    results = []

    for date in dates:
        print(f"Looking for target opportunities on {date}")
        #date_time = timezone.localize(datetime.datetime.combine(date, datetime.time(0, 0)))   

        try:
            astropy_datetime = Time(datetime.datetime.combine(date, datetime.time(12, 0)))

            observation_window_start = observer.sun_set_time(astropy_datetime, which='next', horizon=-12 * u.deg)
            observation_window_end = observer.sun_rise_time(astropy_datetime, which='next', horizon=-12 * u.deg)
            
            observable_minutes = (observation_window_end - observation_window_start).to(u.minute)

            date_is_observable = True

            if observation_window_start.mask or observation_window_end.mask:
                print("No astronomical twilight for this date. Night not valid for observations")
                date_is_observable = False
            elif observable_minutes < min_minutes:
                print(f"{observable_minutes:.0f} obsvervable minute for this night, but {min_minutes:.0f} is required. Night not valid for observations")
                date_is_observable = False

            if date_is_observable:
                observation_window_start_local = observation_window_start.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M')
                observation_window_end_local = observation_window_end.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M')

                print(f"Observing window: {observation_window_start_local} to {observation_window_end_local}")
                print(f"Duration: {observable_minutes:.2f} of astronomical darkness.")

                time_grid = time_grid_from_range([observation_window_start, observation_window_end], time_resolution=10 * u.minute)

                for target in targets:
                    # Evaluate all constraints at once
                    constraint_masks = [
                        constraint(observer, [target], times=time_grid)  # One array of booleans per constraint
                        for constraint in constraints
                    ]

                    # Combine all constraint masks with logical AND
                    combined_mask = constraint_masks[0]
                    for mask in constraint_masks[1:]:
                        combined_mask &= mask

                    target_observable_minutes = combined_mask.sum() * 10 * u.minute

                    
                    if target_observable_minutes >= min_minutes:
                        print(f"Target [{target.name}] is observable on this night")

                        # Find the first time where all constraints are True
                        try:
                            first_index = next(i for i, ok in enumerate(combined_mask) if ok)
                            observation_start_time = time_grid[first_index].to_datetime(timezone=observer.timezone)
                            observation_start_str = observation_start_time.strftime('%Y-%m-%d %H:%M')
                        except StopIteration:
                            observation_start_str = "Unknown"

                        results.append({
                            'date': str(date),
                            'target': target.name,
                            'observable_hours': observable_minutes.to(u.hour).value,
                            'observation_start': observation_start_str
                        })

        except Exception as e:
            print(e)
            continue  # skip problematic days


    if config["make_ics"]:
        print("Making ICS calendar file")
        make_ics(results)
    
    if config['make_calendar']:
        print('Making Calendar Image files')
        make_calendar_images(results, config['calendar_year'])

if __name__ == '__main__':
    main()
