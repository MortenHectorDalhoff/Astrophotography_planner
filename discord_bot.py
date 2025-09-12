from lib.Observer import AstroObserver
from lib.AstroTarget import AstroTarget
import discord
from discord.ext import commands
from datetime import datetime, timedelta
import os
import sys  
import csv

# Define bot intents
intents = discord.Intents.default()

# Initialize the bot
bot = commands.Bot(command_prefix="!", intents=intents)

@bot.event
async def on_ready():

    # Check of log folder is present, if not create it
    if not os.path.exists("logs"):
        os.makedirs("logs")

    await bot.tree.sync()
    print(f"Logged in as {bot.user}")

# Define astro_target_window command
@bot.tree.command(name="astro_target_window", description="Get observation details for a location and target")
async def astro_target_window(
    interaction: discord.Interaction,
    location_name: str,
    target_name: str,
    minimum_observation_minutes: str = "10",
    minimum_altitude: str = "0.0",
    maximum_altitude: str = "90.0",
    minimum_moon_separation: str = "0.0",
    maximum_moon_illumination: str = "100"
    ):

    await interaction.response.defer()

    # start a log file for this interaction
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    username = interaction.user.name.replace(" ", "_")
    log_filename = f"logs/astro_target_window_{timestamp}_{username}.log"
    log_file = open(log_filename, "a", encoding="utf-8")
    console = sys.stdout  # Save the current stdout to restore later
    sys.stdout = log_file

    try:

        # Logging interaction details

        print(f"Received command from {interaction.user}")
        print(f"Location: {location_name}")
        print(f"Target: {target_name}")
        print(f"Minimum Observation Minutes: {minimum_observation_minutes}")
        print(f"Minimum Altitude: {minimum_altitude}")
        print(f"Maximum Altitude: {maximum_altitude}")
        print(f"Minimum Moon Separation: {minimum_moon_separation}")
        print(f"Maximum Moon Illumination: {maximum_moon_illumination}")

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
            minimum_altitude = float(minimum_altitude)
            maximum_altitude = float(maximum_altitude)

            values = (minimum_altitude, maximum_altitude)
            target.add_constraint('altitude', values)
            print(f"Altitude constraints set to min: {minimum_altitude}, max: {maximum_altitude}")
        except ValueError:
            await interaction.response.send_message("Invalid minimum altitude constraint. Use a single value (e.g., 30).")
            return
     
        # add moon separation and illumination constraints
        try:
            values = (float(minimum_moon_separation),)
            target.add_constraint('moon_separation', values)
            print(f"Moon separation constraint set to {minimum_moon_separation} degrees")
        except ValueError:
            await interaction.response.send_message("Invalid moon separation constraint. Use a single value (e.g., 15).")
            return

        
        try:
            values = (float(maximum_moon_illumination),)
            target.add_constraint('moon_illumination', values)
            print(f"Moon illumination constraint set to {maximum_moon_illumination}")
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
            date_str = date.strftime('%Y-%m-%d')
            print(f"Calculating observation window for {date}")

            observation_window = target.get_observation_window(observer, date)
            if observation_window['is_observable']:
                start_str = observation_window['start_str']
                end_str = observation_window['end_str']
                minutes = observation_window['observable_time_str']
                lines.append(f"{date_str:<11} {start_str:<15} {end_str:<15} {minutes:<7}")
            else:
                
                lines.append(f"{date_str:<11} {'-':<15} {'-':<15} {'0 min':<7}")

        table = "```\n" + "\n".join(lines) + "\n```"   

        user_mention = f"<@{interaction.user.id}>"

        full_response = (f"{user_mention}\n\n"
                        f"Observer Location: {observer.location_name}\n"
                        f"Target: {target.pretty_name}\n"
                        f"Minimum Observation Time: {minimum_observation_minutes} minutes\n"
                        f"Altitude: {minimum_altitude} to {maximum_altitude} degrees\n"
                        f"Minimum Moon Separation: {minimum_moon_separation} degrees\n"
                        f"Maximum Moon Illumination: {int(float(maximum_moon_illumination))}%\n\n"
                        f"Observation Windows:\n{table}")

        print("Sending response")
        print(full_response)
        await interaction.followup.send(full_response)

    
    except Exception as e:
        await interaction.followup.send(f"Error: {str(e)}")
        print(f"Error: {str(e)}")

    finally:
        # Restore original stdout and close log file
        sys.stdout = console
        log_file.close()
        print(f"Log saved to {log_filename}")

# Define astro_target_window help command
@bot.tree.command(name="astro_target_window_help", description="Get help information about the bot")
async def help_command(interaction: discord.Interaction):
    help_text = (
        "astro_target_window Bot Commands:\n"
        "/astro_target_window location_name:<location> target_name:<target> [minimum_observation_minutes:<minutes>] "
        "[minimum_altitude:<degrees>] [maximum_altitude:<degrees>] [minimum_moon_separation:<degrees>] [maximum_moon_illumination:<value>]\n\n"
        "Parameters:\n"
        "- location_name: Name of the observer's location (e.g., 'New York, USA')\n"
        "- target_name: Name of the astronomical target (e.g., 'M31')\n"
        "- minimum_observation_minutes: Minimum required observation time in minutes (default: 120)\n"
        "- minimum_altitude: Minimum altitude constraint in degrees (default: 0.0)\n"
        "- maximum_altitude: Maximum altitude constraint in degrees (default: 90.0)\n"
        "- minimum_moon_separation: Minimum angular separation from the Moon in degrees (default: 0.0)\n"
        "- maximum_moon_illumination: Maximum Moon illumination percennt (0 to 100, default: 100)\n\n"
        "Example:\n"
        "/astro_target_window location_name:'New York, USA' target_name:'M31' minimum_observation_minutes:'60' "
        "minimum_altitude:'30' maximum_altitude:'80' minimum_moon_separation:'15' maximum_moon_illumination:'0.5'\n\n"
        "Use /help to see this message again."
    )
    await interaction.response.send_message(help_text)

#Define best_tonight command
@bot.tree.command(name="best_tonight", description="Get the best target to observe tonight from a location")
async def best_tonight(
    interaction: discord.Interaction,
    location_name: str
    ):
    
    await interaction.response.defer()

    # start a log file for this interaction
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    username = interaction.user.name.replace(" ", "_")
    log_filename = f"logs/best_tonight_{timestamp}_{username}.log"
    log_file = open(log_filename, "a", encoding="utf-8")
    console = sys.stdout  # Save the current stdout to restore later
    sys.stdout = log_file

    try:
        # Logginh interaction details

        print(f"Received command from {interaction.user}")
        print(f"Location: {location_name}")

        # Create observer object
        observer = AstroObserver()
        observer.get_location(location_name=location_name)
        print(f"Observer created: {observer}")

        # load targets from csv file into list
        # target,type
        with open('lib/popular_targets.csv', 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            target_list = list(reader)
        

        observerable_targets = []

        for target in target_list:
            target_name = target[0]
            target_type = target[1]

            print(f"Processing target: {target_name} ({target_type})")

            # create target object
            astro_target = AstroTarget()
            astro_target.resolve_target_from_name(name=target_name)
            print(f"Target created: {astro_target.pretty_name}")

            # Minimum 20 dec altitude
            astro_target.add_constraint('altitude', (20.0, 90.0))

            # Minimum 45 dec moon separation
            astro_target.add_constraint('moon_separation', (45.0,))

            # Minimum 60 minutes observation time
            astro_target.set_minimum_observation_minutes(60)

            # Calculate observation window for tonight
            tonight = datetime.now().date() 
            observation_window = astro_target.get_observation_window(observer, tonight)

            if observation_window['is_observable']:
                observerable_targets.append((astro_target.pretty_name, target_type, observation_window['observable_time_str'], observation_window['start_str'], observation_window['end_str']))
                print(f"Target {astro_target.pretty_name} is observable for {observation_window['observable_time_str']}")
        
        # Sort observable targets by observation time
        observerable_targets.sort(key=lambda x: x[2], reverse=True)

        if len(observerable_targets) == 0:
            await interaction.followup.send(f"No observable targets found for {observer.location_name} tonight.")
            return
        
        # Create response message
        print_Nebular = False
        Nebular_table = [
            "Target Name                      Type        Minutes   Start Time      End Time",
            "------------------------------   ----------  -------   --------------- --------------"
        ]
        
        print_galaxy = False
        galaxy_table = [
            "Target Name                      Type        Minutes   Start Time      End Time",
            "------------------------------   ----------  -------   --------------- --------------"
        ]
        
        print_cluster = False
        cluster_table = [
            "Target Name                      Type        Minutes   Start Time      End Time",
            "------------------------------   ----------  -------   --------------- --------------"
        ]
        
        # Show top 10 targets
        for target in observerable_targets:
            target_name = target[0]
            target_type = target[1]
            targert_minutes = target[2]
            target_start = target[3]
            target_end = target[4]

            if target_type == 'Nebular':
                print_Nebular = True
                Nebular_table.append(f"{target_name:<30} {target_type:<10} {targert_minutes:<9} {target_start:<15} {target_end:<15}")

            elif target_type == 'Galaxy':
                print_galaxy = True
                galaxy_table.append(f"{target_name:<30} {target_type:<10} {targert_minutes:<9} {target_start:<15} {target_end:<15}")

            elif target_type == 'Cluster':
                print_cluster = True
                cluster_table.append(f"{target_name:<30} {target_type:<10} {targert_minutes:<9} {target_start:<15} {target_end:<15}")
        
        user_mention = f"<@{interaction.user.id}>"

        if print_Nebular:
            Nebular_message = f"{user_mention}\n```\n" + "\n".join(Nebular_table) + "\n```"
            await interaction.followup.send(Nebular_message)
        if print_galaxy:
            galaxy_message = f"{user_mention}\n```\n" + "\n".join(galaxy_table) + "\n```"
            await interaction.followup.send(galaxy_message)
        if print_cluster:
            cluster_message = f"{user_mention}\n```\n" + "\n".join(cluster_table) + "\n```"
            await interaction.followup.send(cluster_message)

    except Exception as e:
        await interaction.followup.send(f"Error: {str(e)}")
        print(f"Error: {str(e)}")

    finally:
        # Restore original stdout and close log file
        sys.stdout = console
        log_file.close()
        print(f"Log saved to {log_filename}")
        

# Define best_tonight help command
@bot.tree.command(name="best_tonight_help", description="Get help information about the best_tonight command")
async def best_tonight_help_command(interaction: discord.Interaction):
    help_text = (
        "best_tonight Bot Command:\n"
        "/best_tonight location_name:<location>\n\n"
        "Parameters:\n"
        "- location_name: Name of the observer's location (e.g., 'New York, USA')\n\n"
        "Example:\n"
        "/best_tonight location_name:'New York, USA'\n\n"
        "Use /best_tonight_help to see this message again."
    )
    await interaction.response.send_message(help_text)


# Run the bot
if __name__ == "__main__":

    # get the bot token from environment variable
    token = os.environ["DiscordBotToken"]

    # run the bot
    bot.run(token)