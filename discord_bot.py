from lib.Observer import AstroObserver
from lib.AstroTarget import AstroTarget
import discord
from discord.ext import commands
from datetime import datetime, timedelta

# Define bot intents
intents = discord.Intents.default()

# Initialize the bot
bot = commands.Bot(command_prefix="!", intents=intents)

# Define the slash command
@bot.slash_command(name="astroplan", description="Get observation details for a location and target")
async def astroplan(ctx,
                    location_name: str,
                    target_name: str,
                    minimum_observation_minutes: int = 120,
                    min_altitude: float = 0.0,
                    max_altitude: float = 90.0,
                    moon_separation: str = 0.0,
                    moon_illumination: str = 1.0
                    ):
    try:
        # Create observer and target objects
        observer = AstroObserver()
        observer.get_location(location_name=location_name)

        target = AstroTarget()
        target.resolve_target_from_name(name=target_name)

        # Add minimum observation time constraint
        if minimum_observation_minutes:
            try:
                values = int(minimum_observation_minutes)
                target.set_minimum_observation_minutes(values)
            except ValueError:
                await ctx.respond("Invalid minimum observation time constraint. Use a single value in minutes (e.g., 30).")
                return

        # Add altitude constraints
        if min_altitude:
            try:
                min_altitude = float(min_altitude)
            except ValueError:
                await ctx.respond("Invalid minimum altitude constraint. Use a single value (e.g., 30).")
                return


        if max_altitude:
            try:
                max_altitude = float(max_altitude)
            except ValueError:
                await ctx.respond("Invalid maximum altitude constraint. Use a single value (e.g., 90).")
                return

        values = (min_altitude, max_altitude)
        target.add_constraint('altitude', values)

        # add moon separation and illumination constraints
        if moon_separation:
            try:
                values = (float(moon_separation),)
                target.add_constraint('moon_separation', values)
            except ValueError:
                await ctx.respond("Invalid moon separation constraint. Use a single value (e.g., 15).")
                return

        if moon_illumination:
            try:
                values = (float(moon_illumination),)
                target.add_constraint('moon_illumination', values)
            except ValueError:
                await ctx.respond("Invalid moon illumination constraint. Use a single value (e.g., 0.5).")
                return

        # Calculate observation window
        today = datetime.now().date()
        observation_results = []

        for day in range(30):
            date = today + timedelta(days=day)
            observation_window = target.get_observation_window(observer, date)

            if observation_window['is_observable']:
                observation_results.append({
                    'date': date,
                    'start': observation_window['start'],
                    'end': observation_window['end'],
                    'observable_minutes': observation_window['observable_minutes']
                })

        # Prepare Markdown table
        if observation_results:
            table = "| Date       | Start Time       | End Time         | Observable Minutes |\n"
            table += "|------------|------------------|------------------|--------------------|\n"
            for result in observation_results:
                table += f"| {result['date']} | {result['start']} | {result['end']} | {result['observable_minutes']} |\n"
        else:
            table = "No observable windows found in the next 30 days."

        await ctx.respond(f"Observer Location: {observer.location_name}\n"
                          f"Target: {target.pretty_name}\n\n"
                          f"Observation Windows:\n{table}")
    except Exception as e:
        await ctx.respond(f"Error: {str(e)}")

# Run the bot
bot.run("YOUR_BOT_TOKEN")