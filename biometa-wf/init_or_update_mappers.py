"""Initialize or update biological metadata mappers"""
from typing import List
from pathlib import Path
from collections import defaultdict

import yaml
from pymongo import MongoClient


def main():
    attribute_mapper = yaml.safe_load(open("./config/attribute_type.yaml"))
    for attr_class, attr_types in attribute_mapper.items():

        if attr_class == "ignore":
            continue

        class_mapper = get_class_mapper(attr_class)
        for attr_type in attr_types:
            update_class_mapper(attr_type, class_mapper)

        save_class_mapper(attr_class, class_mapper)


def get_class_mapper(attribute_class) -> dict:
    file_name = f"./config/{attribute_class}.yaml"
    if Path(file_name).exists():
        class_mapper = yaml.full_load(open(file_name))
        if hasattr(class_mapper, "get"):
            dd_class_mapper = defaultdict(dict)
            dd_class_mapper.update(class_mapper)
            return dd_class_mapper

    return defaultdict(dict)


def update_class_mapper(attribute_type: str, class_mapper: dict):
    already_processed = class_mapper.get(attribute_type, {}).keys()
    for attr_value in get_attribute_values(attribute_type):
        if attr_value in already_processed:
            continue
        else:
            class_mapper[attribute_type][attr_value] = ""


def get_attribute_values(attribute_type) -> List[str]:
    """Returns a sorted list of attributes give a type."""
    con = db.aggregate(
        [
            {"$unwind": {"path": "$sample.attributes"}},
            {"$match": {"sample.attributes.name": attribute_type}},
            {"$project": {"_id": False, "attr": "$sample.attributes.value"}},
        ]
    )
    return sorted(list(set([x["attr"] for x in con])), key=lambda x: x.lower())


def save_class_mapper(attribute_class: str, class_mapper: defaultdict):
    file_name = f"./config/{attribute_class}.yaml"
    with open(file_name, "w") as fh:
        yaml.dump(dict(class_mapper), fh, allow_unicode=True)


if __name__ == "__main__":
    try:
        client = MongoClient()
        db = client["sramongo"]["ncbi"]
        main()
    finally:
        client.close()
