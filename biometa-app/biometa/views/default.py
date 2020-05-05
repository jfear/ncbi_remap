import json
from collections import defaultdict

from pyramid.view import view_config
from pyramid.request import Request
import pyramid.httpexceptions

from biometa.data_access import (
    get_bioprojects,
    get_bioproject,
    sql_update_biosample,
    get_fly_anatomy,
    get_fly_cell_line,
    get_fly_development,
)

SEX = ["male", "female", "mixed", "None"]
TISSUE = get_fly_anatomy()
DEV_STAGE = get_fly_development()
CELL_TYPE = get_fly_cell_line()
BIOPROJECTS = [x["_id"] for x in get_bioprojects()]


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
        "bioprojects": get_bioprojects(limit=limit, skip=start),
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
            "dev_values": DEV_STAGE,
            "tissue_values": TISSUE,
            "cell_values": CELL_TYPE,
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
            metadata["dev"],
            metadata["tissue"],
            metadata["cell"],
            str(int(metadata.get("perturb") == "on")),
            str(int(metadata.get("complete") == "on")),
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
    curr_idx = BIOPROJECTS.index(accn)
    prev_idx = max(0, curr_idx - 1)

    prev_accn = BIOPROJECTS[prev_idx]
    return pyramid.httpexceptions.HTTPFound(location=f"/project/{prev_accn}")
