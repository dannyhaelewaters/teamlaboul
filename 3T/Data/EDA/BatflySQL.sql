-- create the database and tables for the tweet_py application

CREATE DATABASE IF NOT EXISTS Batfly;

-- use the database
USE Batfly;

-- showing columns
SHOW COLUMNS FROM BatflyData;

-- viewing all infected bat flies
SELECT `unique 3T identifier` FROM BatflyData
WHERE `Laboul infection?` = 'infected';

-- viewing bat family distinct
SELECT DISTINCT `Bat species` FROM BatflyData;

-- viewing infected batflies from laboul infection? column
SELECT DISTINCT Country, `Bat species`, `Laboul infection?`
FROM BatflyData
WHERE `Laboul infection?` = 'infected';

-- obtaining valuecounts of each batfly species that infected with laboul
SELECT `Bat species`, COUNT(`Laboul infection?`) AS `Laboul infection count`
FROM BatflyData
WHERE `Laboul infection?` = 'infected'
GROUP BY `Bat species`
ORDER BY `Laboul infection count` DESC;

-- Obtaining countries of each batfly species that are infected
SELECT Country, `Bat species`, COUNT(`Laboul infection?`) AS `Laboul infection count`
FROM BatflyData
WHERE `Laboul infection?` = 'infected'
GROUP BY `Bat species`, Country
ORDER BY `Laboul infection count` DESC;

-- how many infected batflies within a country
SELECT Country, COUNT(`Laboul infection?`) AS `Laboul infection count`
FROM BatflyData
WHERE `Laboul infection?` = 'infected'
GROUP BY  Country
ORDER BY `Laboul infection count` DESC;


-- How many infected batflies were male or female
SELECT `Sex`, COUNT(`Laboul infection?`) AS `Laboul infection count`
FROM BatflyData
WHERE `Laboul infection?` = 'infected'
GROUP BY  `Sex`
ORDER BY `Laboul infection count` DESC;

-- updating Male as M for Sex column
UPDATE BatflyData
SET `Sex` = 'M'
WHERE
    `Sex` = 'Male';



