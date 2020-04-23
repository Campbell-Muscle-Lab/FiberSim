---
title: JSON Dump Structure
parent: Overview
grand_parent: Core Model
nav_order: 3
---

# The JSON Dump Structure
{:.no_toc}

* TOC
{:toc}

## What Exactly is a JSON?

FiberSim uses a file structure called JSON as it's primary means of model input/ouput. JSON is an easy-to-use data-interchange format that is simple and human-readable. This is much preferred over traditional `txt` files since JSON files enforce a structure that leads them to be easier for both human and machine to interpret.

## The JSON Structure

The JSON format consists of two basic structures:
  + A collection of name and value pairs where the name is used to index the value. For example, the name `"jared"` could index the value `7`. If you're familiar with Python, this is similar to a dictionary.
  + An ordered list of values. This is much like an array, list, or vector in many programming languages.

These basic structures take the forms of:

  + An *object* is an unordered set of the first structure of name and value pairs. Each object begins with a left brace, `{`, and ends with a right brace, `}`. Every name in the object is followed by a colon, then the value. Each name and value pair is seperated by a comma. An example of an object has been included below:
      ```json
      {
        "Jared": 7,
        "Elise": 1200,
        "Allen": 34,
        "Borris": "strings can also be values!",
        "Agatha": null,
        "Harold": true
      }
      ```
  + An *array* is an ordered collection of values. Arrays begin with a left bracket, `[`, and ends with a right bracket, `]`. Values in arrays are separated by commas. The following is an example of an array inside of an object:
      ```json
      {
        "my_array": [10, 11, 12, 13, 14]
      }
      ```
  + A *value* can be a *string* in double quotes; a *number*, meaning an integer or floating point value; a *boolean*, meaning `true` or `false`; or an *object or an *array*. These structures can be nested, meaning you can have arrays of arrays of objects and so forth.

## Where are JSON Files Used in FiberSim?

We use JSON files for the simulation options file, the model file, half-sarcomere status files, and for the visualization instruction files.