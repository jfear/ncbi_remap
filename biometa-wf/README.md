# Biological Meatadata Workflow

The goal of this workflow is to normalize biological metadata to make it more useful. Unfortunately, this requires a lot of iterative hand processing and is not an automated workflow.

## 1. Map Attribute Types (a.k.a. classes/columns)

I want to focus on sex, tissue, cell type, and developmental stage. Unfortunately, users have uploaded hundreds of different feature types with multiple types mapping to these four categories. The first step is to run `./attribute_types.ipynb` and hand map these different columns into `./config/attribute_type.yaml`. The last cell in this notebook copies an attribute type to the clipboard and shows a set of values that are found in this column. You then paste the attribute type into the correct region of the yaml file and it will automatically move to the next attribute.
