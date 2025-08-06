import warnings
import re

from astroplan import FixedTarget
from astroplan.utils import time_grid_from_range
from astroplan.constraints import Constraint
from astroplan.constraints import AltitudeConstraint, AtNightConstraint, MoonSeparationConstraint, MoonIlluminationConstraint

from astropy.coordinates import SkyCoord
from astropy.coordinates import NonRotationTransformationWarning
import astropy.units as u

from astroquery.simbad import Simbad

class AstroTarget:
    def __init__(self):
        self.observable_dates = []
        self.minimum_observation_minutes = 0
        self.constraints = [AtNightConstraint.twilight_astronomical()]

    def __repr__(self):
        """
        Return a string representation of the AstroTarget object.

        Returns
        -------
        str
            A string representation of the AstroTarget object, including its name, RA, and Dec.
        """
        # Return a string representation of the AstroTarget object
        return f"AstroTarget(name={self.name}, ra={self.ra_str}, dec={self.dec_str})"

    def add_constraint(self, constraint_type, values: tuple = None):
        """
        Add a constraint to the target based on the type and values provided.

        Parameters
        ----------
        :param constraint_type: str
            The type of constraint to add. Supported types are 'altitude', 'moon_separation', and 'moon_illumination'.

        :param values: tuple
            A tuple of values required for the constraint. The expected values depend on the constraint type:
            - For 'altitude': one or two values (min, max) in degrees.
            - For 'moon_separation': one value (min) in degrees.
            - For 'moon_illumination': one value (max) as a fraction (0 to 1).

        Raises
        ------
        ValueError
            If the constraint type is unknown or if the number of values provided does not match the expected count
            for the specified constraint type.

        Notes
        -----
        This method allows you to add constraints to the AstroTarget object, which can be used later to check
        whether the target is observable under the specified conditions. The constraints are stored in the
        `self.constraints` list, which can be used with astroplan's observation methods.

        Examples
        --------
        >>> target = AstroTarget('M31')
        >>> target.add_constraint('altitude', (30, 90))  # Set altitude constraint between 30 and 90 degrees
        >>> target.add_constraint('moon_separation', (15,))  # Set moon separation constraint with a minimum of 15 degrees
        >>> target.add_constraint('moon_illumination', (0.5,))  # Set moon illumination constraint with a maximum of 50%
        """
        # Add a constraint to the target based on the type and values provided.
        match constraint_type:

            case 'altitude':
                # Altitude constraint requires one or two values (min, max)
                if len(values) == 2:
                    self.constraints.append(AltitudeConstraint(min=values[0] * u.deg, max=values[1] * u.deg))
                elif len(values) == 1:
                    self.constraints.append(AltitudeConstraint(min=values[0] * u.deg))
                else:
                    raise ValueError("Altitude constraint requires either one or two values (min, max)")

            case 'moon_separation':
                # Moon separation constraint requires one value (min)
                if len(values) == 1:
                    self.constraints.append(MoonSeparationConstraint(min=values[0] * u.deg))
                else:
                    raise ValueError("Moon separation constraint requires one value (min)")

            case 'moon_illumination':
                # Moon illumination constraint requires one value (max)
                if len(values) == 1:
                    self.constraints.append(MoonIlluminationConstraint(max=values[0]/100))
                else:
                    raise ValueError("Moon illumination constraint requires one value (max)")

            case _:
                raise ValueError(f"Unknown constraint type: {constraint_type}")

    def clear_constraints(self):
        """
        Clear all constraints from the target.

        This method removes all constraints that have been added to the target, allowing you to start fresh
        with a new set of constraints if needed.

        Notes
        -----
        This is useful when you want to reset the constraints for the target without creating a new instance.

        Examples
        --------
        >>> target = AstroTarget('M31')
        >>> target.add_constraint('altitude', (30, 90))
        >>> target.clear_constraints()
        """

        # Clear all constraints from the target.
        self.constraints = [AtNightConstraint.twilight_astronomical()]

    def set_minimum_observation_minutes(self, min_time):
        """
        Set the minimum observation time for this target in minutes.

        Parameters
        ----------
        min_time : int
            The minimum observation time in minutes required for this target.

        Raises
        ------
        ValueError
            If the provided minimum time is not a positive integer.

        Notes
        This method allows you to specify the minimum amount of time that must be spent observing this target
        for it to be considered observable. This is useful for filtering out targets that do not meet
        the required observation time criteria.

        Examples
        --------
        >>> target = AstroTarget('M31')
        >>> target.set_minimum_observation_minutes(30)  # Set minimum observation time to 30 minutes
        """
        self.minimum_observation_minutes = min_time

    def get_observation_window(self, observer, date):
        """
        Check if the target is observable from the observer's location on a given date.

        Parameters
        ----------
        observer : Observer
            An instance of the Observer class representing the observer's location and time zone.

        date : datetime.date
            The date for which to check the observability of the target.

        Returns
        -------
        dict
            A dictionary containing the start and end times of the observation window, whether the target is observable,
            and the total observable time in minutes.

        Raises
        ------
        ValueError
            If the observer is not an instance of the Observer class or if the date is not a valid datetime.date object.
        """

        # Get the observation window for the observer on a given date.
        observation_window = observer.get_observation_window(date)

        # Extract start and end times from the observation window
        if observation_window['is_observable'] is False:
            return {
                'start': None,
                'end': None,
                'is_observable': False,
                'observable_minutes': 0 * u.minute
            }

        # Calcualte a list of times where the targer is checked for observability
        time_grid = time_grid_from_range([observation_window_start, observation_window_end],
                                         time_resolution=10 * u.minute)

        # Check if the target is observable at each time in the grid
        constraint_masks = [
            constraint(observer, [self.astroplan_target], times=time_grid)  # One array of booleans per constraint
            for constraint in self.constraints
        ]
        # Combine the masks from all constraints
        combined_mask = constraint_masks[0]
        for mask in constraint_masks[1:]:
            combined_mask &= mask

        # Calculate the total observable time in minutes
        target_observable_minutes = combined_mask.sum() * 10 * u.minute

        if target_observable_minutes >= self.minimum_observation_minutes:

            is_observable = True
            # Find the first time where all constraints are True
            try:
                first_index = next(i for i, ok in enumerate(combined_mask) if ok)
                observation_window_start = time_grid[first_index]
            except StopIteration:
                observation_window_start = "Unknown"

            # Find the last time where all constraints are True
            try:
                last_index = next(i for i in range(len(combined_mask) - 1, -1, -1) if combined_mask[i])
                observation_window_end = time_grid[last_index]
            except StopIteration:
                observation_window_end = "Unknown"
        else:
            is_observable = False


        return_dict = {}
        return_dict['start'] = observation_window_start
        return_dict['end'] = observation_window_end
        return_dict['is_observable'] = is_observable
        return_dict['observable_minutes'] = target_observable_minutes

        return return_dict

    def __get_pretty_name_from_ids(self, ids_string, main_id):
        """
        Extract a human-readable name from the Simbad IDs string.

        Parameters
        ----------
        ids_string : str
            A string containing Simbad IDs, typically in the format 'NAME <pretty_name> | ...'.

        main_id : str
            The main identifier of the target, which is used to filter out the pretty name.

        Returns
        -------
        str
            A human-readable name for the target, derived from the Simbad IDs string.

        Raises
        ------
        ValueError
            If the IDs string is empty or does not contain a pretty name.
        """

        # Extract pretty name from Simbad IDs string.
        aliases = [alias.strip() for alias in ids_string.split('|')]

        # Remove 'NAME ' prefix from aliases
        name_aliases = [
            alias.replace('NAME ', '')
            for alias in aliases
            if alias.startswith('NAME ') and main_id.replace(' ', '') not in alias.replace(' ', '')
        ]
        # Remove 'NAME ' prefix from main_id
        if name_aliases:
            return name_aliases[0]

        return main_id  # No pretty name found

    def resolve_target_from_name(self, name):
        """
        Resolve the target name using Simbad and return the main identifier.

       Parameters
       ----------
       name : str
           The name of the astronomical target, e.g. 'M31' or 'Andromeda Galaxy'.

       Attributes
       ----------
       name : str
           The name of the target.

       observable_dates : list
           List of dates when the target is observable.

       minimum_observation_minutes : int
           Minimum observation time in minutes required for this target.

       constraints : list
           List of constraints for observing the target, such as altitude, moon separation, etc.

       astroplan_target : FixedTarget
           An astroplan FixedTarget object representing the target.

       ra : astropy.units.Quantity
           Right Ascension of the target in degrees.

       dec : astropy.units.Quantity
           Declination of the target in degrees.

       coordinate : SkyCoord
           Astropy SkyCoord object representing the target's coordinates.

       ra_str : str
           Right Ascension formatted as a string in hour angle format.

       dec_str : str
           Declination formatted as a string in degrees, arcminutes, and arcseconds.

       pretty_name : str
           A more human-readable name for the target, derived from Simbad IDs.

       Raises
       ------
       ValueError
           If the target cannot be resolved in Simbad.

       Notes
       -----
       This class is designed to represent an astronomical target and provide methods to check its observability
       from a given observer's location. It uses the Simbad database to resolve target names and retrieve their
       coordinates. The class also allows adding constraints for observing the target, such as altitude and moon separation.

       """

        # Initial variables
        self.name = name

        # Query Simbad for target information
        warnings.filterwarnings('ignore', category=NonRotationTransformationWarning)

        # Add necessary fields to Simbad query
        Simbad.add_votable_fields('ra(d)', 'dec(d)', 'ids')

        # Query Simbad for the target
        result = Simbad.query_object(self.name)
        if result is None:
            raise ValueError(f"Could not resolve target: {self.name}")

        # Extract main_id from the result
        main_id = result['main_id'][0].decode('utf-8') if isinstance(result['main_id'][0], bytes) else \
            result['main_id'][0]

        # Clean up main_id deleting multiple white spaces and 'NAME ' prefix
        main_id = re.sub(r'\s+', ' ', main_id).strip().replace('NAME ', '')  # Clean multiple white spaces

        ids = result['ids'][0]
        self.pretty_name = self.__get_pretty_name_from_ids(ids, main_id)

        # Create SkyCoord object
        self.ra = result['ra'][0] * u.deg
        self.dec = result['dec'][0] * u.deg
        self.coordinate = SkyCoord(ra=self.ra, dec=self.dec)

        # Create FixedTarget with main_id as name
        self.astroplan_target = FixedTarget(name=self.name, coord=self.coordinate)

        # Format RA
        self.ra_str = self.astroplan_target.coord.ra.to_string(unit='hourangle', sep=':')

        # Format Dec with degree, arcmin, arcsec
        deg, arcmin, arcsec = self.astroplan_target.coord.dec.dms
        sign = '+' if deg >= 0 else '-'
        self.dec_str = f"{sign}{abs(deg):.0f}° {abs(arcmin):.0f}' {abs(arcsec):.2f}\""














