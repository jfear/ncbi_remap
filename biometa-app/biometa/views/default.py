from collections import defaultdict

from pyramid.view import view_config
from pyramid.request import Request
import pyramid.httpexceptions

from biometa.data_access import (
    get_bioprojects,
    get_bioproject,
    query_term,
    sql_update_biosample,
    get_fly_anatomy,
    get_fly_cell_line,
    get_fly_development,
    get_biosamples_in_sql,
)

IGNORE_SAMPLE_IN_DB = True  # Set if you want to focus on samples not yet annotated
IGNORE_SAMPLES = []
if IGNORE_SAMPLE_IN_DB:
    IGNORE_SAMPLES = get_biosamples_in_sql()

SEX = ["male", "female", "mixed", "None"]
TISSUE = get_fly_anatomy()
DEVEL_STAGE = get_fly_development()
CELL_LINE = get_fly_cell_line()
BIOPROJECTS = [x["_id"] for x in get_bioprojects(ignore_already_processed=IGNORE_SAMPLES)]
print(len(BIOPROJECTS))


# /
@view_config(route_name="home")
def home(request):
    return pyramid.httpexceptions.HTTPFound(location="/0")


# /{start}
@view_config(route_name="bioprojects", renderer="../templates/bioprojects.pt")
def bioprojects_page(request):
    # Simple pagination
    limit = 20
    start = int(request.matchdict.get("start"))
    _next = int(start) + limit
    _prev = max(0, int(start) - limit)

    return {
        "bioprojects": get_bioprojects(
            limit=limit, skip=start, ignore_already_processed=IGNORE_SAMPLES
        ),
        "next": str(_next),
        "prev": str(_prev),
    }


# /project/{accn}
@view_config(
    route_name="bioproject", renderer="../templates/bioproject.pt", request_method="GET",
)
def bioproject_page(request: Request):
    accn = request.matchdict.get("accn")
    bioproject = get_bioproject(accn)

    if bioproject:
        return {
            "bioproject": bioproject,
            "sex_values": SEX,
            "devel_values": DEVEL_STAGE,
            "tissue_values": TISSUE,
            "cell_values": CELL_LINE,
        }

    raise pyramid.httpexceptions.HTTPNotFound()


@view_config(
    route_name="bioproject", renderer="../templates/bioproject.pt", request_method="POST",
)
def bioproject_post(request: Request):
    # Process from
    result_dict = defaultdict(dict)
    for k, v in request.POST.items():
        biosample, metadata_type = k.strip().split("_")
        result_dict[biosample][metadata_type] = v

    # Prep for SQL UPSERT
    results = [
        (
            biosample,
            metadata["sex"],
            metadata["devel"],
            metadata["tissue"],
            metadata["cell"],
            metadata["notes"],
            str(int(metadata.get("genetic") == "on")),
            str(int(metadata.get("diet") == "on")),
            str(int(metadata.get("chemical") == "on")),
            str(int(metadata.get("radiation") == "on")),
            str(int(metadata.get("temperature") == "on")),
            str(int(metadata.get("other") == "on")),
            str(int(metadata.get("control") == "on")),
        )
        for biosample, metadata in result_dict.items()
    ]

    sql_update_biosample(results)

    # Reload bioproject
    curr_accn = request.matchdict.get("accn")
    return pyramid.httpexceptions.HTTPFound(location=f"/project/{curr_accn}")


# /next/{accn}
@view_config(route_name="next_project")
def next_project(request: Request):
    accn = request.matchdict.get("accn")
    curr_idx = BIOPROJECTS.index(accn)
    next_idx = curr_idx + 1

    if next_idx > len(BIOPROJECTS):
        next_idx = 0

    next_accn = BIOPROJECTS[next_idx]
    return pyramid.httpexceptions.HTTPFound(location=f"/project/{next_accn}")


# /prev/{accn}
@view_config(route_name="previous_project")
def previous_project(request: Request):
    accn = request.matchdict.get("accn")
    curr_idx = BIOPROJECTS.index(accn) or 0
    prev_idx = max(0, curr_idx - 1)

    prev_accn = BIOPROJECTS[prev_idx]
    return pyramid.httpexceptions.HTTPFound(location=f"/project/{prev_accn}")


# /search/{col}/{term}
@view_config(
    route_name="query", renderer="../templates/bioproject.pt", request_method="GET",
)
def query_page(request: Request):
    col = request.matchdict.get("col")
    term = request.matchdict.get("term")

    bioproject = query_term(col, term)

    if bioproject:
        return {
            "bioproject": bioproject,
            "sex_values": SEX,
            "devel_values": DEVEL_STAGE,
            "tissue_values": TISSUE,
            "cell_values": CELL_LINE,
        }

    raise pyramid.httpexceptions.HTTPNotFound()


@view_config(
    route_name="query", renderer="../templates/bioproject.pt", request_method="POST",
)
def query_post(request: Request):
    # Process from
    result_dict = defaultdict(dict)
    for k, v in request.POST.items():
        biosample, metadata_type = k.strip().split("_")
        result_dict[biosample][metadata_type] = v

    # Prep for SQL UPSERT
    results = [
        (
            biosample,
            metadata["sex"],
            metadata["devel"],
            metadata["tissue"],
            metadata["cell"],
            metadata["notes"],
            str(int(metadata.get("genetic") == "on")),
            str(int(metadata.get("diet") == "on")),
            str(int(metadata.get("chemical") == "on")),
            str(int(metadata.get("radiation") == "on")),
            str(int(metadata.get("temperature") == "on")),
            str(int(metadata.get("other") == "on")),
            str(int(metadata.get("control") == "on")),
        )
        for biosample, metadata in result_dict.items()
    ]

    sql_update_biosample(results)

    # Reload bioproject
    curr_col = request.matchdict.get("col")
    curr_term = request.matchdict.get("term")
    return pyramid.httpexceptions.HTTPFound(location=f"/search/{curr_col}/{curr_term}")
