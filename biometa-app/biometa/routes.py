def includeme(config):
    config.add_static_view("static", "static", cache_max_age=3600)
    config.add_route("home", "/")
    config.add_route("next_project", "/next/{accn}")
    config.add_route("previous_project", "/prev/{accn}")
    config.add_route("bioproject", "/project/{accn}")
    config.add_route("bioprojects", "/{start}")
