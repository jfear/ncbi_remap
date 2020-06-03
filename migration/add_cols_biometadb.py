# %%
import sqlite3

# %%
con = sqlite3.connect("../data/biometa_old.db")
cur = con.cursor()
cur.execute("SELECT * FROM biometa")
dat = cur.fetchall()
con.close()
# %%
transformed_dat = []
for row in dat:
    transformed_dat.append((
        row[0], # biosample
        row[1], # sex
        row[2], # devel_stage
        row[3], # tissue
        row[4], # cell_line
        "", # notes
        "0", # genetic
        "0", # diet
        "0", # chemical
        "0", # radiation
        "0", # temperature
        str(row[6]), # other
        "0", # control
    ))


# %%
transformed_dat[0]

# %%
con = sqlite3.connect("../data/biometa.db")
cur = con.cursor()
sql_upsert = """INSERT INTO biometa(biosample, sex, devel_stage, tissue, cell_line, notes, genetic, diet, chemical, radiation, temperature, other, control)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ON CONFLICT(biosample) DO UPDATE SET
        sex = excluded.sex,
        devel_stage = excluded.devel_stage,
        tissue = excluded.tissue,
        cell_line = excluded.cell_line,
        notes = excluded.notes,
        genetic = excluded.genetic,
        diet = excluded.diet,
        chemical = excluded.chemical,
        radiation = excluded.radiation,
        temperature = excluded.temperature,
        other = excluded.other,
        control = excluded.control
"""
cur.executemany(sql_upsert, transformed_dat)
con.commit()
con.close()

# %%
