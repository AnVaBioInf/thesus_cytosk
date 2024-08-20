library("recount3")

# downloading available projects
human_projects <- recount3::available_projects()

# subsetting only sra projects with more than 10 samples
sra=subset(human_projects, file_source == "sra" & n_samples>=10 & project_type == "data_sources")[,'project']
print(sra)

# writing sra ids to file
write(sra, file = "sra.txt", sep = " ")

#—-------------------
  
# #bash
# # installing bio package for sta metadata parsing 
# pip3 install bio
# 
# # for each sra id extracting metadata, and than choosing only lines where there is a word develop*, and if so, printing sra id as well, all the output writing to a file
# data=$(cat sra.txt)
# for variable in $data; do bio search $variable -all | grep -m1 "develop*" && echo $variable; done > develop_all.txt
# for variable in $data; do bio search $variable -all | grep -m1 "breast*" && echo $variable; done > breast_all.txt
