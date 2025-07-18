import datetime
from geopy.geocoders import Nominatim
from timezonefinder import TimezoneFinder
from astroplan import Observer
import astropy.units as u
from astropy.time import Time

class Observer:

    def __init__(self, location_name):
        """
        Initialize an observer at a given location by name.

        Parameters
        ----------
        location_name : str
            The name of the location (e.g., "New York, USA").

        Attributes:
        ----------
        name : str
            The name of the observer.

        latitude : float
            The latitude of the observer's location in degrees.

        longitude : float
            The longitude of the observer's location in degrees.

        timezone : str
            The timezone of the observer's location.

        Raises
        ------
        ValueError
            If the location cannot be found or if latitude/longitude cannot be determined.

        Example
        -------
        >>> observer = Observer("New York, USA")
        >>> print(observer.name)
        Ovserver(name=New York, USA, latitude=40.7128, longitude=-74.0060, timezone=America/New_York)
        """

        # Get lat and lon from name
        loc = Nominatim(user_agent="GetLoc")
        geo_location = loc.geocode(location_name)
        self.latitude = geo_location.latitude
        self.longitude = geo_location.longitude

        # Look up timezone from lat adn lon
        tz_finder = TimezoneFinder()
        self.timezone = tz_finder.timezone_at(lat=latitude, lng=longitude)

    def __repr__(self):
        """
        Return a string representation of the observer.

        Returns
        -------
        str
            A string representation of the observer including name, latitude, longitude, and timezone.
        """
        return f"Ovserver(name={self.location_name}, latitude={self.latitude}, longitude={self.longitude}, timezone={self.timezone})"

    def get_observer(self):
        """
        Get the observer object for the observer's location.

        Returns
        -------
        Observer
            An astroplan Observer object initialized with the observer's latitude, longitude, name, and timezone.
        """
        # Create an observer object using astroplan
        return Observer(latitude=self.latitude, longitude=self.longitude, name=self.location_name, timezone=self.timezone)

    def get_observation_window(self, date):
        """
        Get the observation window for the observer on a given date.

        Parameters
        ----------
        date : datetime.date
            The date for which to calculate the observation window.

        Returns
        -------
        dict
            A dictionary containing the start and end times of the observation window, whether it is observable,
            and the number of observable minutes.

        Raises
        ------
        ValueError
            If the date is not a valid datetime.date object.
        """

        observer = self.get_observer()
        astropy_datetime = Time(datetime.datetime.combine(date, datetime.time(12, 0)))

        # Calculate the start and end of astronomical twilight
        observation_window_start = observer.sun_set_time(astropy_datetime, which='next', horizon=-12 * u.deg)
        observation_window_end = observer.sun_rise_time(astropy_datetime, which='next', horizon=-12 * u.deg)

        # Calculate the observable minutes
        observable_minutes = (observation_window_end - observation_window_start).to(u.minute)

        # Check if the observation window is valid
        if observation_window_start.mask or observation_window_end.mask:
            print("No astronomical twilight for this date. Night not valid for observations")
            is_observable = False
        else:
            is_observable = True

        return_dict = {}
        return_dict['start'] = observation_window_start
        return_dict['end'] = observation_window_end
        return_dict['is_observable'] = is_observable
        return_dict['observable_minutes'] = observable_minutes

        return return_dict
