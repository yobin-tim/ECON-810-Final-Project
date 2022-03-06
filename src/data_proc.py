from matplotlib.pyplot import axis
import pandas as pd
import numpy as np

def read_data(path):
    '''
    Reads in the data from the given path or url.
    '''
    data = pd.read_stata(path)
    return data

def filter_data(data, start_year, end_year, variables, filters ):
    '''
    Filters the data for the given start and end year.
    Inputs:
        data: pandas dataframe
        start_year: num
        end_year: num
        variables: list or dict
        filters: dict {variable: filter} if len(filters) == 1 filter using '==' if len(filters) == 2 filter using '<=' and '>=' 
    Returns:
        df: pandas dataframe
    '''
    # First we drop individuals from the SEO oversample
    df = data[data.x11104LL == 'Main Sample    11']
    # Then we filter the dataframe by by and years
    df = df[(df.year >= start_year) & (df.year <= end_year)]
    
    if isinstance(variables, list):
        # If the variables are a list, we return the dataframe with the variables
        df = df[variables]
    elif isinstance(variables, dict):
        # If the variables are a dict, we return the dataframe with the variables and rename the columns
        # Create a new dataframe with the data we want
        df = df[variables.keys()]
        # Rename the columns
        df.rename(columns=variables, inplace=True)
    
    for variable, filter in filters.items():
        # For each variable and filter, we filter the dataframe
        if len(filter) == 1:
            # If the filter is a single number, we filter the dataframe with '=='
            df = df[df[variable] == filter[0]]
        elif len(filter) == 2:
            # If the filter is a list of two numbers, we filter the dataframe with '<=' and '>='
            df = df[(df[variable] >= filter[0]) & (df[variable] <= filter[1])]
import pandas as pd
import numpy as np

def read_data(path):import pandas as pd
import numpy as np

def read_data(path):
    '''
    Reads in the data from the given path or url.
    '''
    data = pd.read_stata(path)
    return data

def filter_data(data, start_year, end_year, variables, filters ):
    '''
    Filters the data for the given start and end year.
    Inputs:
        data: pandas dataframe
        start_year: num
        end_year: num
        variables: list or dict
        filters: dict {variable: filter} if len(filters) == 1 filter using '==' if len(filters) == 2 filter using '<=' and '>=' 
    Returns:
        df: pandas dataframe
    '''
    # First we drop individuals from the SEO oversample
    df = data[data.x11104LL == 'Main Sample    11']
    # Then we filter the dataframe by by and years
    df = df[(df.year >= start_year) & (df.year <= end_year)]
    
    if isinstance(variables, list):
        # If the variables are a list, we return the dataframe with the variables
        df = df[variables]
    elif isinstance(variables, dict):
        # If the variables are a dict, we return the dataframe with the variables and rename the columns
        # Create a new dataframe with the data we want
        df = df[variables.keys()]
        # Rename the columns
        df.rename(columns=variables, inplace=True)
    
    for variable, filter in filters.items():
        # For each variable and filter, we filter the dataframe
        if len(filter) == 1:
            # If the filter is a single number, we filter the dataframe with '=='
            df = df[df[variable] == filter[0]]
        elif len(filter) == 2:
            # If the filter is a list of two numbers, we filter the dataframe with '<=' and '>='
            df = df[(df[variable] >= filter[0]) & (df[variable] <= filter[1])]
import pandas as pd
import numpy as np

def read_data(path):
    '''
    Reads in the data from the given path or url.
    '''
    data = pd.read_stata(path)
    return data

def filter_data(data, start_year, end_year, variables, filters ):
    '''
    Filters the data for the given start and end year.
    Inputs:
        data: pandas dataframe
        start_year: num
        end_year: num
        variables: list or dict
        filters: dict {variable: filter} if len(filters) == 1 filter using '==' if len(filters) == 2 filter using '<=' and '>=' 
    Returns:
        df: pandas dataframe
    '''
    # First we drop individuals from the SEO oversample
    df = data[data.x11104LL == 'Main Sample    11']
    # Then we filter the dataframe by by and years
    df = df[(df.year >= start_year) & (df.year <= end_year)]
    
    if isinstance(variables, list):
        # If the variables are a list, we return the dataframe with the variables
        df = df[variables]
    elif isinstance(variables, dict):
        # If the variables are a dict, we return the dataframe with the variables and rename the columns
        # Create a new dataframe with the data we want
        df = df[variables.keys()]
        # Rename the columns
        df.rename(columns=variables, inplace=True)
    
    for variable, filter in filters.items():
        # For each variable and filter, we filter the dataframe
        if len(filter) == 1:
            # If the filter is a single number, we filter the dataframe with '=='
            df = df[df[variable] == filter[0]]
        elif len(filter) == 2:
            # If the filter is a list of two numbers, we filter the dataframe with '<=' and '>='
            df = df[(df[variable] >= filter[0]) & (df[variable] <= filter[1])]

    return df 

def filter_by_earnigs(data, earnings_variable, earnings_cuttof, number_of_years):
    '''
    Filters the data for the given earnings_variable and earnings_cuttof and number_of_years.
    Inputs:
        data: pandas dataframe
        earnings_variable: str
        earnings_cuttof: num or list if num then  earnings_cuttof = [earnings_cuttof, inf]
        number_of_years: num
    Outputs:
        df: pandas dataframe
    '''

    if not isinstance(earnings_cuttof, list):
        earnings_cuttof = [earnings_cuttof, np.inf]

    # Next we create a summary for any individual how many years has earned more than the lower bound cutoff
    summary = data.groupby("id").agg({earnings_variable : lambda column : column[(column >= earnings_cuttof[0])].count()})
    # Then we select the individuals who fall between inside the range for the number of years
    id_s = summary[summary[earnings_variable] >= number_of_years].index
    # Set Personal id as the index
    df = data.set_index("id")
    # Then we filter the dataframe by the summary
    df = df.loc[id_s]

    # Repeat for the upper bound cutoff but remove individual hwo have earned more than the upper bound cutoff at least one year
    summary = df.groupby("id").agg({earnings_variable : lambda column : column[(column <= earnings_cuttof[1])].count()})
    id_s = summary[summary[earnings_variable] >= 1].index
    # df = df.loc[id_s]
    print(id_s)
    return df

    return df 

def filter_by_earnigs(data, earnings_variable, earnings_cuttof, number_of_years):
    '''
    Filters the data for the given earnings_variable and earnings_cuttof and number_of_years.
    Inputs:
        data: pandas dataframe
        earnings_variable: str
        earnings_cuttof: num or list if num then  earnings_cuttof = [earnings_cuttof, inf]
        number_of_years: num
    Outputs:
        df: pandas dataframe
    '''

    if not isinstance(earnings_cuttof, list):
        earnings_cuttof = [earnings_cuttof, np.inf]

    # Next we create a summary for any individual how many years has earned more than the lower bound cutoff
    summary = data.groupby("id").agg({earnings_variable : lambda column : column[(column >= earnings_cuttof[0])].count()})
    # Then we select the individuals who fall between inside the range for the number of years
    id_s = summary[summary[earnings_variable] >= number_of_years].index
    # Set Personal id as the index
    df = data.set_index("id")
    # Then we filter the dataframe by the summary
    df = df.loc[id_s]

    # Repeat for the upper bound cutoff but remove individual hwo have earned more than the upper bound cutoff at least one year
    summary = df.groupby("id").agg({earnings_variable : lambda column : column[(column <= earnings_cuttof[1])].count()})
    id_s = summary[summary[earnings_variable] >= 1].index
    # df = df.loc[id_s]
    print(id_s)
    return df

    '''
    Reads in the data from the given path or url.
    '''
    data = pd.read_stata(path)
    return data

def filter_data(data, start_year, end_year, variables, filters ):
    '''
    Filters the data for the given start and end year.
    Inputs:
        data: pandas dataframe
        start_year: num
        end_year: num
        variables: list or dict
        filters: dict {variable: filter} if len(filters) == 1 filter using '==' if len(filters) == 2 filter using '<=' and '>=' 
    Returns:
        df: pandas dataframe
    '''
    # First we drop individuals from the SEO oversample
    df = data[data.x11104LL == 'Main Sample    11']
    # Then we filter the dataframe by by and years
    df = df[(df.year >= start_year) & (df.year <= end_year)]
    
    if isinstance(variables, list):
        # If the variables are a list, we return the dataframe with the variables
        df = df[variables]
    elif isinstance(variables, dict):
        # If the variables are a dict, we return the dataframe with the variables and rename the columns
        # Create a new dataframe with the data we want
        df = df[variables.keys()]
        # Rename the columns
        df.rename(columns=variables, inplace=True)
    
    for variable, filter in filters.items():
        # For each variable and filter, we filter the dataframe
        if len(filter) == 1:
            # If the filter is a single number, we filter the dataframe with '=='
            df = df[df[variable] == filter[0]]
        elif len(filter) == 2:
            # If the filter is a list of two numbers, we filter the dataframe with '<=' and '>='
            df = df[(df[variable] >= filter[0]) & (df[variable] <= filter[1])]

    return df 

def filter_by_earnigs(data, earnings_variable, earnings_cuttof, number_of_years):
    '''
    Filters the data for the given earnings_variable and earnings_cuttof and number_of_years.
    Inputs:
        data: pandas dataframe
        earnings_variable: str
        earnings_cuttof: num or list if num then  earnings_cuttof = [earnings_cuttof, inf]
        number_of_years: num
    Outputs:
        df: pandas dataframe
    '''

    if not isinstance(earnings_cuttof, list):
        earnings_cuttof = [earnings_cuttof, np.inf]

    # Next we create a summary for any individual how many years has earned more than the lower bound cutoff
    summary = data.groupby("id").agg({earnings_variable : lambda column : column[(column >= earnings_cuttof[0])].count()})
    # Then we select the individuals who fall between inside the range for the number of years
    id_s = summary[summary[earnings_variable] >= number_of_years].index
    # Set Personal id as the index
    df = data.set_index("id")
    # Then we filter the dataframe by the summary
    df = df.loc[id_s]

    # Repeat for the upper bound cutoff but remove individual hwo have earned more than the upper bound cutoff at least one year
    summary = df.groupby("id").agg({earnings_variable : lambda column : column[(column <= earnings_cuttof[1])].count()})
    id_s = summary[summary[earnings_variable] >= 1].index
    # df = df.loc[id_s]
    print(id_s)
    return df

    return df 

def filter_by_earnigs(data, earnings_variable, earnings_cuttof, number_of_years):
    '''
    Filters the data for the given earnings_variable and earnings_cuttof and number_of_years.
    Inputs:
        data: pandas dataframe
        earnings_variable: str
        earnings_cuttof: num or list if num then  earnings_cuttof = [earnings_cuttof, inf]
        number_of_years: num
    Outputs:
        df: pandas dataframe
    '''

    if not isinstance(earnings_cuttof, list):
        earnings_cuttof = [earnings_cuttof, np.inf]

    # Next we create a summary for any individual how many years has earned more than the lower bound cutoff
    summary = data.groupby("id").agg({earnings_variable : lambda column : column[(column >= earnings_cuttof[0])].count()})
    # Then we select the individuals who fall between inside the range for the number of years
    id_s = summary[summary[earnings_variable] >= number_of_years].index
    # Set Personal id as the index
    df = data.set_index("id")
    # Then we filter the dataframe by the summary
    df = df.loc[id_s]

    # Repeat for the upper bound cutoff but remove individual hwo have earned more than the upper bound cutoff at least one year
    summary = df.groupby("id").agg({earnings_variable : lambda column : column[(column > earnings_cuttof[1])].count()})
    id_s = summary[summary[earnings_variable] >= 1].index
    df = df.drop(id_s, axis=0)

    # Finally remove individuals who have negative earnings on any year
    summary = df.groupby("id").agg({earnings_variable : lambda column : column[(column < 0 )].count()})
    id_s = summary[summary[earnings_variable] >= 1].index
    df = df.drop(id_s, axis=0)

    return df
