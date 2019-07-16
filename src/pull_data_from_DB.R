# This file uses dbplyr package to fetch filtered and grouped tables from the DB.
# Here client side operations should be keeped to bare minimum and no server side
# operations. All the new files created in the workspace should be displayed.
# These variables should be saved so the can be releoaded faster. A query_db
# boolean variable can be used to switch between DB fetch requests and reading
# from an RDS for faster data reload.
