Automated-Image-Analysis-for-Ki67-Scoring-in-Gastroenteropancreatic-Neuroendocrine-Tumors

The ML and non-ML analyses each created a file containing the coordinates of their corresponding detected cells. After the manual count (ground truth), we also saved a simple file with the coordinates of manually detected cells. Here, the code reads the numbers and coordinates of total and positive tumor cells detected by both analysis tools and those from the ground truth. Then, it performs a distance calculation between the tumor cells detected by each analysis tool and the ground truth to find the best match. The purpose is to identify the correct tumor cell detections. At the end, a Word file is generated containing all the results.


