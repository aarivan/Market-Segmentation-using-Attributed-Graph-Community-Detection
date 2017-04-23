# Market-Segmentation-using-Attributed-Graph-Community-Detection

# Overview: 
Market segmentation divides a broad target market into subsets of consumers or businesses that have or are perceived to have common needs, interests, and priorities. These segments help firms or businesses focus on their target groups effectively and allocate resources efficiently. Traditional segmentation methods are solely based on attribute data such as demographics (age, sex, ethnicity, education, etc.) and psychographic profiles (lifestyle, personality, motives, etc.). However, social networks have recently become important for marketing. Depending on the nature of the market, social relations can even become vital in forming segments. Such social relations combined with demographic properties can be used to find more relevant subsets of consumers or businesses (i.e., communities)

In this project, we aim to find such market segments given social network data. These social relations can be captured in a graph framework where nodes represent customers/users and edges represent some social relationship. The properties belonging to each customer/user can be treated as node attributes. Hence, market segmentation becomes the problem of community detection over attributed graphs, where the communities are formed based on graph structure as well as attribute similarities.

# Influence Propagation: 
It has been discovered that by targeting the influential users/groups, desirable marketing goals can be achieved. Therefore, one way to evaluate the quality of the market segments (communities) is to influence an entity in each segment and measure how fast the influence propagates over the entire network. The faster that influence propagates through the entire network, the more likely an advertising campaign, for example, will be successful.

# Goal: 
To implement a community detection algorithm for attributed graphs, find the relevant market segments, and evaluate the obtained segments via influence propagation.

# Datasets:
We are provided with two small datasets: fb_caltech_small_edgelist.txt and fb_caltech_small_attrlist.csv. 

This dataset contains a facebook network of a US university (given as an edgelist) with each node corresponding to a user profile having the following attributes: student/faculty status, gender, major, second major, dorm, and year information. For the similarity convenience, these attribute values have been converted into asymmetric binary variables.

# Project Requirements
We are expected to implement an algorithm from the paper “Community detection based on structural and attribute similarities.”

# Cluster Evaluation
The code will measure and compare the time steps taken by the influence propagation algorithm to achieve maximum propagation and compare this with kmeans clustering.

Note: This script will assumes you have a file named communities.txt
