## Adding a new parameter to CIAlign

If you wish to add a new parameter to CIAlign, please use the following steps - or let us know if you have any problems with any of them so we can help.

To add a parameter:

* If the parameter is a float or int, add a name, default, minimum and maximum value to columns 1, 2, 3 and 4 of CIAlign/ranges.txt
* Add the parameter to the CIAlign/argP.py help documentation, following the style of the other parameters
* Add the parameter the function call in CIAlign/runCIAlign.py
* Edit the function itself, including adding the parameter to the function definition and describing it in the docstring, add comments if needed to explain your steps
* Edit the unit tests for the appropriate module file to include the new parameter
* Run unit tests in the root directoy as `python -m unittest tests/*py` and ensure that they pass

To edit the manual:
* Edit man/user_guide_template.md, following the style for the other parameters and adding your parameter. The minimum, maximum and default will be added automatically if you use this style.
* Run `python man/update_manual.py`, which should update 
man/user_guide.md, man/user_guide.pdf and README.md
* Commit your changes and make a pull request
