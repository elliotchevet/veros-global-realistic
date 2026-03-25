from veros.settings import Setting

SETTINGS = dict(
    "start_date": Setting("",str, "Start date of simulation"),
    "dt_forcing": Setting(3, int, "Time step for forcing in hours"),
    )
