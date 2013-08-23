-- MySQL dump 10.9
--
-- Host: localhost    Database: ciwogGCAG
-- ------------------------------------------------------
-- Server version	4.1.12

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `ClusterIntron2Gene`
--

DROP TABLE IF EXISTS `ClusterIntron2Gene`;
CREATE TABLE `ClusterIntron2Gene` (
  `intronNum` int(10) default NULL,
  `annotationUID` varchar(100) NOT NULL default '',
  `clusterUID` varchar(20) NOT NULL default '',
  `clusterIntronUID` int(10) NOT NULL default '0',
  `phase` char(1) default NULL,
  `intronType` varchar(10) default NULL,
  `org` varchar(25) NOT NULL default '',
  `alnPos` varchar(10) default NULL,
  `proteinPos` varchar(10) default NULL,
  `intronUID` int(10) default NULL,
  `altgroup` int(10) default NULL,
  `mod_intronType` varchar(10) default NULL,
  `termDn` varchar(4) default NULL,
  `intronSeq` text,
  `start` int(50),
  `stop` int(50),
  `segment` int(50),
  `strand` char(1),
  KEY clusterUID (clusterUID),
  KEY clusterIntronUID (clusterIntronUID),
  KEY org (org),
  KEY annotationUID (annotationUID),
  KEY `intronType` (`intronType`),
  KEY `intronUID` (`intronUID`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--
-- Table structure for table `ClusterIntrons`
--

DROP TABLE IF EXISTS `ClusterIntrons`;
CREATE TABLE `ClusterIntrons` (
  `clusterIntronUID` int(10) NOT NULL default '0',
  `clusterUID` varchar(20) NOT NULL default '',
  `alnPosStart` char(10) default NULL,
  `alnPosEnd` char(10) default NULL,
  PRIMARY KEY  (`clusterIntronUID`,`clusterUID`),
  KEY `alnPosStart` (`alnPosStart`),
  KEY `alnPosEnd` (`alnPosEnd`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--
-- Table structure for table `Clusters`
--

DROP TABLE IF EXISTS `Clusters`;
CREATE TABLE `Clusters` (
  `clusterUID` varchar(20) NOT NULL default '',
  `distinctIntronCount` int(10) default NULL,
  `geneCount` int(10) default NULL,
  `avgIntronAmount` decimal(8,4) default '-9999.9999',
  `intronRange` text,
  `orgCount` int(10) default NULL,
  `distinctOrgL` text,
  `u12GeneCount` int(10) default NULL,
  `u12MultipleGeneCount` int(10) default NULL,
  `avgSim` decimal(8,4) default '-9999.9999',
  `avgCov` decimal(8,4) default '-9999.9999',
  `minSim` decimal(8,4) default '-9999.9999',
  `maxSim` decimal(8,4) default '-9999.9999',
  `maxCov` decimal(8,4) default '-9999.9999',
  `minCov` decimal(8,4) default '-9999.9999',
  `u12OrgL` text,
  `u12IntronCount` int(10) default NULL,
  `u12WdualIntronCount` int(10) default NULL,
  `u12WdualOrgL` text,
  `u12WdualGeneCount` int(10) default NULL,
  `dualIntronCount` int(10) default NULL,
  `dualGeneCount` int(10) default NULL,
  `u12OrgCount` int(10) default NULL,
  `dualOrgCount` int(10) default NULL,
  `u12WdualOrgCount` int(10) default NULL,
  `dualOrgL` text,
  `u12hGeneCount` int(10) default NULL,
  `u12hOrgL` text,
  `u12hIntronCount` int(10) default NULL,
  `u12hOrgCount` int(10) default NULL,
  PRIMARY KEY  (`clusterUID`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

--
-- Table structure for table `clusterMembers`
--

DROP TABLE IF EXISTS `clusterMembers`;
CREATE TABLE `clusterMembers` (
  `clusterUID` varchar(20) NOT NULL default '',
  `org` varchar(6) default NULL,
  `annotationUID` int(10) default NULL,
  `geneName` varchar(45) default NULL,
  `altgroup` int(10) default NULL,
  `needReview` tinyint(1) default NULL,
  `redundant` tinyint(1) default NULL,
  `intronLess` tinyint(1) default '0',
  UNIQUE KEY `org_2` (`org`,`annotationUID`),
  KEY `org` (`org`),
  KEY `annotationUID` (`annotationUID`),
  KEY `altgroup` (`altgroup`),
  KEY `clusterUID` (`clusterUID`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

