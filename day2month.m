function month_name = day2month(day)
% Converts day of the year to month
% Input: day - day of the year (integer)
% Output: month - corresponding month (integer)

% Define the number of days in each month (non-leap year)
days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% Check if input is valid
if day < 1 || day > 365
    error('Invalid input. Day must be between 1 and 365.');
end

% Initialize variables
total_days = 0;
month = 1;

% Loop through months to find the corresponding month
while total_days < day
    total_days = total_days + days_per_month(month);
    month = month + 1;
end

% Adjust for overshoot
if total_days > day
    month = month - 1;
end

if month == 1
    month_name = 'January';
elseif month == 2
        month_name = 'February';
elseif month == 3
        month_name = 'March';
elseif month == 4
        month_name = 'April';
elseif month == 5
        month_name = 'May';
elseif month == 6
        month_name = 'June';
elseif month == 7
        month_name = 'July';
elseif month == 8
        month_name = 'August';
elseif month == 9
        month_name = 'September';
elseif month == 10
        month_name = 'October';
elseif month == 11
        month_name = 'November';
elseif month == 12
        month_name = 'December';
end

end
