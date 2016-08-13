# Angular-Regression
Run an angular regression model with your dataset.


## How to run the app
You must have R and the Shiny package installed, and then running the app is as easy as entering the following command:
shiny::runGist('https://gist.github.com/aureliennicosia/272fd65592a2fac6422ca7e29eebf0dd')

## Example
Run the app with the dataset 'caribou.txt' with the following specifications:
- the travelled distance is d;
- the direction of the animal is y;
- xcut represents the direction from the animal to the closest regenerating cut;
- xcenter represents the direction from the animal to the closest centroid of previously visited locations;
- the formula is y~xcut+xcenter+xcut:d+xcenter:d .
