from lib.Observer import AstroObserver
from lib.AstroTarget import AstroTarget
import discord
from discord.ext import commands
from datetime import datetime, timedelta
import os

# Define bot intents
intents = discord.Intents.default()

# Initialize the bot
bot = commands.Bot(command_prefix="!", intents=intents)

@bot.event
async def on_ready():
    await bot.tree.sync()
    print(f"Logged in as {bot.user}")

# Define the slash command
@bot.tree.command(name="astroplan", description="Get observation details for a location and target")
async def astroplan(
    interaction: discord.Interaction,
    location_name: str,
    target_name: str,
    minimum_observation_minutes: str = "10",
    min_altitude: str = "0.0",
    max_altitude: str = "90.0",
    moon_separation: str = "0.0",
    moon_illumination: str = "1.0"
    ):

    await interaction.response.defer()

    try:

        # Logginh interaction details

        print(f"Received command from {interaction.user}")
        print(f"Location: {location_name}")
        print(f"Target: {target_name}")
        print(f"Minimum Observation Minutes: {minimum_observation_minutes}")
        print(f"Min Altitude: {min_altitude}")
        print(f"Max Altitude: {max_altitude}")
        print(f"Moon Separation: {moon_separation}")
        print(f"Moon Illumination: {moon_illumination}")

        # Create observer and target objects
        observer = AstroObserver()
        observer.get_location(location_name=location_name)
        print(f"Observer created: {observer}")

        target = AstroTarget()
        target.resolve_target_from_name(name=target_name)
        print(f"Target created: {target.pretty_name}")

        # Add minimum observation time constraint
        try:
            minimum_observation_minutes = int(minimum_observation_minutes)
            target.set_minimum_observation_minutes(minimum_observation_minutes)
            print(f"Minimum observation time constraint set to {minimum_observation_minutes} minutes")
        except ValueError:
            await interaction.response.send_message("Invalid minimum observation time constraint. Use a single value in minutes (e.g., 30).")
            return

        # Add altitude constraints
        try:
            min_altitude = float(min_altitude)
            max_altitude = float(max_altitude)

            values = (min_altitude, max_altitude)
            target.add_constraint('altitude', values)
            print(f"Altitude constraints set to min: {min_altitude}, max: {max_altitude}")
        except ValueError:
            await interaction.response.send_message("Invalid minimum altitude constraint. Use a single value (e.g., 30).")
            return
     
        # add moon separation and illumination constraints
        try:
            values = (float(moon_separation),)
            target.add_constraint('moon_separation', values)
            print(f"Moon separation constraint set to {moon_separation} degrees")
        except ValueError:
            await interaction.response.send_message("Invalid moon separation constraint. Use a single value (e.g., 15).")
            return

        
        try:
            values = (float(moon_illumination),)
            target.add_constraint('moon_illumination', values)
            print(f"Moon illumination constraint set to {moon_illumination}")
        except ValueError:
            await interaction.response.send_message("Invalid moon illumination constraint. Use a single value (e.g., 0.5).")
            return

        # Calculate observation window
        today = datetime.now().date()
        lines = [
                "Date        Start Time      End Time        Minutes",
                "----------- --------------- --------------- -------"
            ]

        for day in range(30):
            date = today + timedelta(days=day)
            print(f"Calculating observation window for {date}")

            observation_window = target.get_observation_window(observer, date)

            if observation_window['is_observable']:
                date_str = date.strftime('%Y-%m-%d')
                start_str = observation_window['start_str']
                end_str = observation_window['end_str']
                minutes = observation_window['observable_time_str']
                lines.append(f"{date_str:<11} {start_str:<15} {end_str:<15} {minutes:<7}")
            else:
                date_str = date.strftime('%Y-%m-%d')
                broken_constraints = observation_window.get('broken_constraints', [])
                lines.append(f"{date_str:<11} {broken_constraints}")

        table = "```\n" + "\n".join(lines) + "\n```"   

        print("Sending response")
        await interaction.followup.send(f"Observer Location: {observer.location_name}\n"
                        f"Target: {target.pretty_name}\n"
                        f"Minimum Observation Time: {minimum_observation_minutes} minutes\n"
                        f"Altitude: {min_altitude} to {max_altitude} degrees\n"
                        f"Minimum Moon Separation: {moon_separation} degrees\n"
                        f"Maximum Moon Illumination: {int(float(moon_illumination)*100)}%\n\n"
                        f"Observation Windows:\n{table}")
    except Exception as e:
        await interaction.response.send_message(f"Error: {str(e)}")

# Define a help command
@bot.tree.command(name="astroplan_help", description="Get help information about the bot")
async def help_command(interaction: discord.Interaction):
    help_text = (
        "AstroPlan Bot Commands:\n"
        "/astroplan location_name:<location> target_name:<target> [minimum_observation_minutes:<minutes>] "
        "[min_altitude:<degrees>] [max_altitude:<degrees>] [moon_separation:<degrees>] [moon_illumination:<value>]\n\n"
        "Parameters:\n"
        "- location_name: Name of the observer's location (e.g., 'New York, USA')\n"
        "- target_name: Name of the astronomical target (e.g., 'M31')\n"
        "- minimum_observation_minutes: Minimum required observation time in minutes (default: 120)\n"
        "- min_altitude: Minimum altitude constraint in degrees (default: 0.0)\n"
        "- max_altitude: Maximum altitude constraint in degrees (default: 90.0)\n"
        "- moon_separation: Minimum angular separation from the Moon in degrees (default: 0.0)\n"
        "- moon_illumination: Maximum Moon illumination fraction (0.0 to 1.0, default: 1.0)\n\n"
        "Example:\n"
        "/astroplan location_name:'New York, USA' target_name:'M31' minimum_observation_minutes:'60' "
        "min_altitude:'30' max_altitude:'80' moon_separation:'15' moon_illumination:'0.5'\n\n"
        "Use /help to see this message again."
    )
    await interaction.response.send_message(help_text)

# Run the bot
if __name__ == "__main__":

    # get the bot token from environment variable
    token = os.environ["DiscordBotToken"]

    # run the bot
    bot.run(token)