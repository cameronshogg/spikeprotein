# STAT 331 - Spike Protein
Spike Protein Accuracy Prediction for COVID-19

The objective of this project was to develop a multiple linear regression model to best predict the accuracy of computer-generated structures for the COVID-19 spike protein. To train the model, a training dataset with 1946 samples was provided. Alongside the response variable, accuracy, there were 685 explanatory variables. In initial analysis of the dataset, two of the explanatory variables were removed due to perfect multicollinearity. An initial test of the effect of removing highly multicollinear variables and observations considered outliers on predicted Root Mean Squared Error (RMSE) led to the removal of 115 more explanatory variables and 206 observations from the dataset. After this, nine models were tested against each other to gain a sense of what model selection techniques would produce the lowest RMSE. The three techniques with the lowest RMSE were tested against each other, trained on 80% of the initial dataset with outliers removed and tested on the remaining 20% of the data in the validation set. The model with the lowest RMSE of 0.571 in this test was selected by an ICM algorithm with a penalty of 2. This same selection technique was used on the entire dataset with the 115 highly multicollinear variables and 206 outliers removed, and the model generated was chosen as the final model used to predict the 1946 accuracy values in the test set.
