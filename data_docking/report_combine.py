#!/share/public/software/Python-2.7.13/bin python2.7
#-*-coding=utf-8-*-
from optparse import OptionParser
import create_family_xml
import create_common_xml
import ConfigParser

import urllib2,urllib

import json  
import copy  
import glob 
import stat
import time
import sys
import os  
import re
global pdf_name
global absolute_json
global ISOTIMEFORMAT
ISOTIMEFORMAT='%Y-%m-%d %X'
absolute_json="./"
pdf_name=""

###统一字段名称
def trans_ne_appendix_name(list):
	result=[]
	if len(list):
		for item in list:
			tmp=item
			tmp['Sample']=item['sample']
			del(tmp['sample'])
			if 'gene' in item.keys():
				tmp['Gene']=item['gene']
				del(tmp['gene'])
			if 'position' in item.keys():
				tmp['Position']=item['position']
				del(tmp['position'])
			if 'NM' in item.keys():
				tmp['Transcript']=item['NM']
				del(tmp['NM'])
			if 'exon' in item.keys():
				tmp['Exon']=item['exon']
				del(tmp['exon'])
			if 'hgvs_c' in item.keys():
				tmp['hgvs.c']=item['hgvs_c']
				del(tmp['hgvs_c'])
			if 'hgvs_p' in item.keys():
				tmp['hgvs.p']=item['hgvs_p']
				del(tmp['hgvs_p'])
			if 'var_type' in item.keys():
				tmp['Type']=item['var_type']
				del(tmp['var_type'])
			##阴性附录的genotype特殊处理
			if 'genotype' in item.keys():
				tmp['Genotype']=item['genotype']
				del(tmp['genotype'])
			###其他家系的杂合性
			if 'father' in item.keys():
				tmp['Father']=item['father']
				del(tmp['father'])
			else:
				tmp['Father']='NA'
			if 'mother' in item.keys():
				tmp['Mother']=item['mother']
				del(tmp['mother'])
			else:
				tmp['Mother']='NA'
			if 'others' in item.keys():
				tmp['Others']=item['others']
				del(tmp['others'])
			else:
				tmp['Others']='NA'
			if 'clinical_level' in item.keys():
				tmp['ACMGLevel']=item['clinical_level']
				del(tmp['clinical_level'])
			if 'disease' in item.keys():
				tmp['Disease']=item['disease']
				del(tmp['disease'])
			if 'inheritance' in item.keys():
				tmp['Inheritance']=item['inheritance']
				del(tmp['inheritance'])
			if 'maf' in item.keys():
				tmp['AlleleFrequency']=item['maf']
				del(tmp['maf'])
			if 'rs' in item.keys():
				tmp['dbSNP']=item['rs']
				del(tmp['rs'])
			if 'REMARK' in item.keys():
				tmp['REMARK']=item['REMARK']
				#del(tmp['rs'])
			result.append(tmp)
	
	return result	
###提取真正的一代验证图片和wescnv图片
def get_verify_imgs_wescnv_imgs(dict,config_path=os.path.split(os.path.realpath(__file__))[0]):
	config = ConfigParser.ConfigParser()
	config_path = os.path.join(config_path, "config.txt")
	config.read(config_path)
	cnv_to_it_prefix= config.get("sample_path", "CNV_TO_IT_PATH")
	####
	wes_verify_imgs={}
	wescnv_imgs_names=''
	num=1
	verify_wescnv_imgs_list=[]
	##1)WES检测范围内的CNV
	tmp_wescnv_imgs_name=[]
	for item in dict.keys():
		tmp=[]
		title=dict[item][0]
		if 'WES检测范围' not in title:
			tmp=dict[item]
			wes_verify_imgs[num]=tmp
			num=num+1
		elif 'WES检测范围' in title:
			image=dict[item][1]
			###拷贝到IT服务器
			command="scp " + image +" "+cnv_to_it_prefix
			os.system(command)
			img_name=dict[item][1].split('/')[-1]
			tmp_wescnv_imgs_name.append(img_name)
			pass
		wescnv_imgs_names=(',').join(tmp_wescnv_imgs_name)
	return wes_verify_imgs,wescnv_imgs_names
##阴性附录疾病描述加粗
def get_bold_appendix(list):
	disease_bold=[]
	for item in list:
		if 'REMARK' in item:
			temp=item['REMARK'].split('。')
			for jtem in temp:
				disease=temp[0].replace('\n','').replace(' ','').split('的临床特征')[0]
				if disease not in disease_bold:
					disease_bold.append(disease)
	return disease_bold
###获取关键词
def get_key_words(sample_code,config_path=os.path.split(os.path.realpath(__file__))[0]):
	result=''
	key_words_list=[]
	config = ConfigParser.ConfigParser()
	config_path = os.path.join(config_path, "config.txt")
	config.read(config_path)
	keywords_prefix= config.get("sample_path", "KEY_WORDS")
	str="wes_website/Phenolyzer/"+sample_code+"/cn_input.file" ##修改文件，2018.12.14
	keywords_path=os.path.join(keywords_prefix,str)
	if os.path.exists(keywords_path):
		f=open(keywords_path)
		lines=f.readlines()
		for line in lines:
			term=line.strip()
			if term not in key_words_list and term !='':
				key_words_list.append(term)
	result='、'.join(key_words_list)
	return result,key_words_list
###获取家系其他成员的核心CNV信息
def get_other_family_members_cnv(sample_id,bus_code_combine,list):
	for id in  bus_code_combine:
		if id == sample_id:
			pass
		else:
			family_info={}
			url=r'http://10.100.11.55/berry/gdwes/getByBusCode/'+id
			try:
				html=urllib2.urlopen(url,timeout=10)
			except  urllib2.URLError,e:
				print "The clinical API did not return data, please check the url:'%s'!"%url
				sys.exit(1)
			sample_info=html.read()
			if sample_info:
				hjson = json.loads(sample_info)	#转为json数据
				result=hjson['success']
				if result==True:
					family_info = hjson
					if len(family_info['vus']):
						for item in family_info['vus'] :
							if 'NA' not in item['chr_position']:
								#print "aaa:",type(item)
								list.append(item)
	return	list
			
#########将NA替换为NA
def trans_NA_to_ND(list):
	result=[]
	##core_data
	for item in list:
		for i in item.keys():
			if 'NA' == item[i]:
				item[i]='ND'
	result=list
	return result
####填充空的扩展信息
def fill_ND_to_Extend(sample_id):
	sites={}
	sites['Gene']="ND"
	sites['Sample']=sample_id
	sites['Type']="ND"
	sites['Position']="ND"
	sites['Transcript']="ND"
	sites['Exon']="ND"
	sites['hgvs.c']="ND"
	sites['hgvs.p']="ND"
	sites['Genotype']="ND"
	sites['dbSNP']="ND"
	sites['AlleleFrequency']="ND"
	sites['Disease']="ND"
	sites['Inheritance']="ND"
	sites['ACMGLevel']="ND"
	sites['Typical_age_of_onset']="ND"
	return sites
####删除散样扩展报告中与核心或提示重复的位点
def del_bulk_repeat_sites(list1,list2,list3,sample_id):
	sites_extend={}
	index=1
	extend_new=[]
	sites_core={}
	sites_note={}
	for item in list1:
		if 'ND' not in item['Gene'] and 'ND' not in item['Position']:
			#gene=item['Gene']
			pos=item['Position']
			tmp=pos
			aa=tmp.encode('unicode-escape').decode('string_escape')
			sites_core[aa]=1
	for jtem in list2:
		if 'ND' not in jtem['Gene'] and 'ND' not in jtem['Position']:
			#gene=jtem['Gene']
			pos=jtem['Position']
			tmp=pos
			bb=tmp.encode('unicode-escape').decode('string_escape')
			sites_note[bb]=1
	tmp_extend=[]
	for item in list3:
		index=0
		if 'ND' not in item['Gene'] and 'ND' not in item['Position']:
			#gene=item['Gene']
			pos=item['Position']
			tmp=pos
			if tmp in sites_core.keys() or tmp in sites_note.keys():
				pass
			else:
				tmp_extend.append(item)
		else:
			tmp_extend.append(item)
	##如果去重后的先证者信息为空了，需要用ND填充
	if len(tmp_extend)==0:
		sites=fill_ND_to_Extend(sample_id)
		extend_new.append(sites)
	else:
		for jtem in tmp_extend:
			extend_new.append(jtem)
	return extend_new


####删除家系扩展报告中与核心或提示重复的位点
def del_family_repeat_sites(list1,list2,list3,sample,relationship,couple_id):
	sample_id=sample
	relation=relationship
	sites_core={}
	sites_extend={}
	index=1
	extend_new=[]
	sites_note={}
	for item in list1:
		if 'ND' not in item['Gene'] and 'ND' not in item['Position']:
			gene=item['Gene']
			pos=item['Position']
			tmp=gene+pos
			aa=tmp.encode('unicode-escape').decode('string_escape')
			sites_core[aa]=1
	for jtem in list2:
		if 'ND' not in jtem['Gene'] and 'ND' not in jtem['Position']:
			gene=jtem['Gene']
			pos=jtem['Position']
			tmp=gene+pos
			bb=tmp.encode('unicode-escape').decode('string_escape')
			sites_note[bb]=1
	
	###1.夫妻特殊处理；2.普通家系，只对先证者去重，且要保证先证者，父亲、母亲按照之前的顺序
	if ('夫' in relation) or ('妻' in relation ):
		extend_new=del_couple_repeat_sites(sites_core,sites_note,list3,sample,couple_id)
	else:
		tmp_extend=[]
		##1.处理先证者
		for item in list3:
			if item['Sample']==sample_id:
				if 'ND' not in item['Gene'] and 'ND' not in item['Position']:
					gene=item['Gene']
					pos=item['Position']
					tmp=pos
					if tmp in sites_core.keys() or tmp in sites_note.keys():
						pass
					else:
						tmp_extend.append(item)
				else:
					tmp_extend.append(item)
			else:
				pass
		##如果去重后的先证者信息为空了，需要用ND填充
		if len(tmp_extend)==0:
			sites=fill_ND_to_Extend(sample_id)
			extend_new.append(sites)
		else:
			for jtem in tmp_extend:
				extend_new.append(jtem)
		##2.加入非先证者
		for item in list3:
			if item['Sample']!=sample_id:
				extend_new.append(item)			
	return extend_new
######去除表格里的多余空行
def del_NA_line(dict):
	tmp=[]
	for item in dict:
		##防止位点信息中包含NA，如NACLP
		if "ND" in item[0] and "ND" in item[1]:
			tmp.append(item)
		for i in tmp:
			dict.remove(i)
	return dict


	
#--------------获取科诺安cnv图片路径信息--------------------_#
def get_cnv_picture_result(cnv_picture_file=[],bus_code="",status="",config_path=os.path.split(os.path.realpath(__file__))[0]): 
	config = ConfigParser.ConfigParser()
	config_path = os.path.join(config_path, "config.txt")
	config.read(config_path)
	cnv_from_it_prefix= config.get("sample_path", "CNV_FROM_IT_PATH")
	cnv_bit_path_prefix=config.get("sample_path", "CNV_BIT_PATH")
	cnv_to_it_prefix= config.get("sample_path", "CNV_TO_IT_PATH")
	#利用sftp获取CNV结果图
	if status.strip():
		#command="sftp -oPort=3033 "+ "wes@10.100.16.45:/zonghe/sharedisk/sharedisk/cnv/newImage/"+bus_code.strip()+"*/*.png  /share/work1/zhanglj/cnv_picture/"
		command="sftp -oPort=3033 "+ cnv_from_it_prefix + bus_code.strip()+"*/*.png "+cnv_bit_path_prefix
		os.system(command)
	#path="/share/work1/zhanglj/cnv_picture/"+bus_code.strip()+"*"
	path=cnv_bit_path_prefix + bus_code.strip()+"*"
	file_list=glob.glob(path)
	if file_list:
		cnv_picture_file=file_list
		cnv_picture_file=file_list
		for item in cnv_picture_file:
			##传给IT
			#command="scp "+item +" gps@10.100.21.50:/zonghe/sharedisk/gps-shared/WES/cnv_picture/"
			command="scp "+item +" " +cnv_to_it_prefix
			print "SCP KNA:",command
			os.system(command)
	else:
		pass
	if not len(cnv_picture_file):
		cnv_picture_file=[]
	return cnv_picture_file 


#返回拼接好之后的质控路径信息
def get_all_qc_path(bus_code_list,config_path=os.path.split(os.path.realpath(__file__))[0]):
	config = ConfigParser.ConfigParser()
	config_path = os.path.join(config_path, "config.txt")
	config.read(config_path)
	qc_prefix= config.get("sample_path", "QC_PATH")
	extend_prefix=config.get("sample_path","EXTEND_PATH")
	qc_path_list=[]
	result=""
	if len(bus_code_list):
		for bus_code in bus_code_list:
			#str="*/analysis/"+bus_code+"/6.QC/QC/"+bus_code+"_QC.xls"
			str="*/analysis/"+bus_code+"/6.QC/QC/"+bus_code+".Core.gene.xls" ##修改文件，2018.12.14
			
			sample_path=os.path.join(qc_prefix,str)
			if len(glob.glob(sample_path))>0:
				qc_path_list.append(glob.glob(sample_path)[0])
			else:
				print "The QC file path: %s is not exists!"%(sample_path)
				sys.exit(1)
		if len(qc_path_list)==1:
			result=qc_path_list[0]
		elif len(qc_path_list) >1:
			result=",".join(qc_path_list)
		else:
			print "There is no qc path for any bus code!"
			sys.exit(1)
	else:
		print "There is no bus code information!"
		sys.exit(1)
	return result
###########获取扩展报告绝对路径
def get_absolute_extend_path(bus_code_list,config_path=os.path.split(os.path.realpath(__file__))[0]):
	config = ConfigParser.ConfigParser()
	config_path = os.path.join(config_path, "config.txt")
	config.read(config_path)
	extend_prefix=config.get("sample_path", "EXTEND_PATH")
	extend_path_list=[]
	result=""
	if len(bus_code_list):
		for bus_code in bus_code_list:
			str_extend="*/analysis/"+bus_code+"/5.Interpretation/ACMG/"+bus_code+"_Extended.xls"
			tmp_extend_path=os.path.join(extend_prefix,str_extend)
			##扩展绝对路径
			if len(glob.glob(tmp_extend_path))>0:
				if "_10G" in glob.glob(tmp_extend_path)[0]:
					extend_path_list.append(glob.glob(tmp_extend_path)[1])
				else:
					extend_path_list.append(glob.glob(tmp_extend_path)[0])
			else:
				print "The Extend file path: %s is not exists!"%(tmp_extend_path)
				sys.exit(1)
		if len(extend_path_list)==1:
			result=extend_path_list[0]
		elif len(extend_path_list) >1:
			result=",".join(extend_path_list)
		else:
			print "There is no extend path for any bus code!"
			sys.exit(1)
	else:
		print "There is no bus code information!"
		sys.exit(1)
	return result
###########获取qc的信息########
def get_all_qc_data(qc_path_string):
	qc_list=qc_path_string.split(",")
	
	result=[]
	for item in qc_list:
		sites={}
		sample_id=item.split("/")[-1].split(".Core.gene.xls")[0]
		sites["sample"] = sample_id
		sites['seq_project']=u'人类全外显子组'
		if os.path.exists(item):
			f=open(item)
			
			lines=f.readlines()
			for line in lines:
				if "PCT_TARGET_BASES_20X" in line:
					tmp_cover_20x=float(line.split()[-1].strip())*100
					tmp=format(tmp_cover_20x,'.2f')
					sites['cover_20x'] = str(tmp)+"%"
		else:
			print "The QC data %s is not exists"%(item)
			sys.exit(1)
		result.append(sites)
	return result
def deal_every_file_extend_site(site_line,sampleid):
	result=[]
	if len(site_line):
		sites={}
		line_data=site_line.split("\t")
		item=[data.strip() for data in line_data]
		sites['Gene']=item[0]
		sites['Sample']=sampleid
		sites['Type']=item[1]
		sites['Position']=item[2]
		sites['Transcript']=item[3]
		sites['Exon']=item[4]
		sites['hgvs.c']=item[5]
		sites['hgvs.p']=item[6]
		if item[7]=='het':
			sites['Genotype']=u'杂合'
		elif item[7]=='hom':
			sites['Genotype']=u'纯合'
		elif item[7]=='hem':
			sites['Genotype']=u'半合子'
		sites['dbSNP']=item[8]
		sites['AlleleFrequency']=item[9]
		sites['Disease']=item[10]
		sites['Inheritance']=item[11]
		sites['ACMGLevel']=item[12]
		sites['Typical_age_of_onset']=item[13]
	return sites
def get_all_extend_data(extend_path_string):
	extend_list=extend_path_string.split(",")
	sites=[]
	for item in extend_list:
		if os.path.exists(item):
			sampleid=item.split("_Extended.xls")[0].split("/")[-1]
			f=open(item)
			lines=f.readlines()
			sample_site={}
			if len(lines) >1:
				for line in lines[1:]:
					sites.append(deal_every_file_extend_site(line,sampleid))
			else:
				sample_site['Sample']=sampleid
				sample_site['Gene']='ND'
				sample_site['Type']='ND'
				sample_site['Position']='ND'
				sample_site['Transcript']='ND'
				sample_site['Exon']='ND'
				sample_site['hgvs.c']='ND'
				sample_site['hgvs.p']='ND'
				sample_site['Genotype']='ND'
				sample_site['dbSNP']='ND'
				sample_site['AlleleFrequency']='ND'
				sample_site['Disease']='ND'
				sample_site['Inheritance']='ND'
				sample_site['ACMGLevel']='ND'
				sample_site['Typical_age_of_onset']='ND'
				sites.append(sample_site)
		else:
			print "The Extend table %s is not exists"%(extend_file)
			sys.exit(1)
	return sites
	
##############获取GD接口的位点信息#############
def get_gd_info(subnumber,outfile):
	sample_code=subnumber.split('-')[0].replace('Z','')
	url=r'http://10.100.11.55/berry/gdwes/getByBusCode/'+sample_code
	try:
		html=urllib2.urlopen(url,timeout=10)
	except  urllib2.URLError,e:
		print "The clinical API did not return data, please check the url:'%s'!"%url
		sys.exit(1)
	sample_info=html.read()
	exam_gd={}
	exam_gd['subnumber']=''
	sample_list=[]
	exam_gd['qc_path']=[]
	exam_gd['qc_data']=[]
	exam_gd['extend_path']=[]
	exam_gd['extend_data']=[]
	exam_gd['articles']=[]
	exam_gd['conclusion']=[]
	exam_gd['verifyResult']=[]
	exam_gd['check_result']=''
	exam_gd['key_words']='' ##关键词，add20190619
	exam_gd['appendix_red']='' #add20190619
	exam_gd['appendix_bold']='' #add20190619
	input_info=""
	temp=""
	
	pattern="[-_a-zA-Z0-9]+"
	bus_code_combine=[]
	if sample_info:
		hjson = json.loads(sample_info)	#转为json数据
		result=hjson['success']
		if result==True:
			exam_gd = hjson
	exam_gd['subnumber']=subnumber
	###暂定使用report，补充报告
	if len(exam_gd['supplement_report']):
		###判断补充类型
		for item in exam_gd['supplement_report']:
			if item['subnumber']!=subnumber:
				exam_gd['supplement_report']=[]
				break
			####判断重采学的情况C0
			if sample_code in item['sample'] :
				if item['sample'].endswith('C0'):
					item['sample']=sample_code
			if item['sample']!=sample_code:
				continue
			if 'supplement_type' in item.keys():
				exam_gd['supplement_reason']=item['supplement_type'].encode('utf-8')
				supplement_type=item['supplement_type'].encode('utf-8')
				if u'补充原检测报告中的核心位点验证' in supplement_type or '对变异位点信息补充验证结果' in supplement_type:
					exam_gd['supplement_type']='add_sup_sites'
				elif u'由于新增临床信息进行再分析' in supplement_type or u'其他家属加做全外进行再分析' in supplement_type:
					exam_gd['supplement_type']='alter_reanalysis_wes'
				exam_gd['original_clinicalinfo']=item['original_clinicalinfo']
				exam_gd['add_clinicalinfo']=item['add_clinicalinfo']
				exam_gd['supplement_conclusion']=item['supplement_conclusion']
	if len(exam_gd['supplement_report']):
		for item in exam_gd['supplement_report']:
			if '\n' in item['family_carry']:
				tmp=item['family_carry'].replace('\n','')
				item['family_carry']=tmp
			else:
				pass
	
	###补充阴性附录的情况，需要排除其他补充类型
	if len(exam_gd['negative_appendix']):
		##判断真实的补充类型，可能是正常报告阴性附录
		for item in exam_gd['negative_appendix']:
			if 'subnumber' not in item.keys():
				continue
			if item['subnumber']!=subnumber:
				exam_gd['negative_appendix']=[]
				break
			####判断重采学的情况C0
			if sample_code in item['sample'] :
				if item['sample'].endswith('C0'):
					item['sample']=sample_code
			if item['sample']!=sample_code:
				continue
			if 'original_clinicalinfo' not in item.keys():
				continue
			if 'supplement_type' in item.keys():
				exam_gd['supplement_reason']=item['supplement_type']
				supplement_type=item['supplement_type'].encode('utf-8')
				if u'增加阴性附录报告' in supplement_type:
					exam_gd['supplement_type']='sup_ne_appendix'
				exam_gd['original_clinicalinfo']=item['original_clinicalinfo']
				exam_gd['add_clinicalinfo']=item['add_clinicalinfo']
				exam_gd['supplement_conclusion']=item['supplement_conclusion']
	if len(exam_gd['negative_appendix']):
		for item in exam_gd['negative_appendix']:
			if '\n' in item['family_carry']:
				tmp=item['family_carry'].replace('\n','')
				item['family_carry']=tmp
			else:
				pass
	print 'SUB_TYPE:', exam_gd['supplement_type']	
	####增加家系成员的顺序，保证报告中的家系信息也是按照这个顺序
	##add20190214
	exam_gd['relation_order']=[u'先证者',u'父亲',u'父親',u'母親',u'母亲',u'夫妻',u'丈夫',u'妻子',u'爷爷',u'奶奶',u'外公',u'外婆',u'女儿',u'儿子',u'哥哥',u'大哥',u'二哥',u'三哥',u'姐姐',u'大姐',u'二姐',u'三姐',u'弟弟',u'大弟',u'二弟',u'三弟',u'妹妹',u'大妹',u'二妹',u'三妹',u'小妹',u'姑母',u'大姑',u'二姑',u'三姑',u'小姑',u'叔叔',u'大叔',u'二叔',u'三叔',u'小叔',u'伯伯',u'大伯',u'二伯',u'三伯',u'舅舅',u'大舅',u'二舅',u'三舅',u'小舅',u'姨',u'大姨',u'二姨',u'三姨',u'小姨',u'堂兄',u'堂弟',u'堂姐',u'堂妹',u'表兄',u'表弟',u'表姐',u'表妹',u'侄子',u'侄女',u'外甥',u'外甥女',u'孙子',u'孙女',u'外孙子',u'外孙女',u'胎儿1',u'胎儿2',u'弟弟或妹妹',u'姑父']
	order_relation_data=[]
	exam_gd['relationship']=""
	if len(hjson['clinic']):
		for item in hjson['clinic']:
			if item['bus_code']==sample_code:
				temp = item['family_test']
				sample=item
	else:
		print "The  Sample:%s is not exists in GD database, please check it!"%(sample_code)
		sys.exit(1)
	
	if temp:
		if "，" in temp:
			temp = temp.replace("，", ",")
		if "," in temp:
			temp_data = temp.split(",")
			####按照家系的顺序存取数据
			for key in exam_gd['relation_order']:
				for item in temp_data:
					##外甥也在外甥女之内，所以不能用 if key in ,须用==
					if key == item.split(' ')[0]:
						order_relation_data.append(item)
			temp_data=order_relation_data
			if len(temp_data):
				temp_data = [item.strip() for item in temp_data]
				#####增加判断，客服现在录入家系时会有空格，如FU 8E1048FU0
				for item in temp_data:
					if ' ' in item:
						new_item=item.split(' ')[1]
						search=re.search(pattern,new_item)
						temp_sample=search.group()
						relation=item.split(' ')[0]
						if("夫" in relation) or ("妻" in relation):
							exam_gd['relationship']=relation
						else:
							exam_gd['relationship']=''
					else:
						search=re.search(pattern,item)
						temp_sample=search.group()
					##去掉末尾带R0的编号，包括FUR0，MUR0,add20190521
					str_end=temp_sample.rfind('R0')
					if str_end!=-1:
						continue
					####
					if ("V" in temp_sample) or ("MV" in temp_sample) or ("XV" in temp_sample) or ("XV1" in temp_sample) or ("XV2" in temp_sample) or ("XV3" in temp_sample) or ('-T' in temp_sample):
						continue
					bus_code_combine.append(temp_sample)
		else:
			search=re.search(pattern,temp)
			one_temp_sample=search.group()
			bus_code_combine=[one_temp_sample]
	##如果sample_code中有C0时，需要把family_test中的原编号忽略add20190521
	for item in bus_code_combine:
		if 'C0' in item and sample_code == item:
			tmp=item.replace('C0','')
			if tmp in bus_code_combine:
				bus_code_combine.remove(tmp)
		elif 'C0' in item:
			tmp=item.replace('C0','0')
			if tmp in bus_code_combine:
				bus_code_combine.remove(tmp)
	exam_gd['bus_code_combine'] = bus_code_combine
	##新增，关键词
	exam_gd['key_words'],exam_gd['appendix_red']=get_key_words(sample_code)
	exam_gd['appendix_bold']=get_bold_appendix(exam_gd['negative_appendix'])
	###只有改做wes+再分析时才需要关注result
	if exam_gd['supplement_type']=='add_sup_sites' or exam_gd['supplement_type']=='alter_reanalysis_wes':
		if hjson['supplement_report'][0]['result']=="negative":
			exam_gd['check_result']=u'阴性'
		elif hjson['supplement_report'][0]['result']=="positive":
			exam_gd['check_result']=u'阳性'
		elif hjson['supplement_report'][0]['result']=="unknown":
			exam_gd['check_result']=u'未知'
		else:
			exam_gd['check_result']=hjson['supplement_report'][0]['result']
	elif exam_gd['supplement_type']=='sup_ne_appendix':
		if hjson['negative_appendix'][0]['result']=="negative":
			exam_gd['check_result']=u'阴性'
		elif hjson['negative_appendix'][0]['result']=="positive":
			exam_gd['check_result']=u'阳性'
		elif hjson['negative_appendix'][0]['result']=="unknown":
			exam_gd['check_result']=u'未知'
		else:
			exam_gd['check_result']=hjson['negative_appendix'][0]['result']
	else:
		print "The  Sample:%s does not has Clear supplement types, please check it!"%(sample_code)
		sys.exit(1)	
	#-----10 WES+CNV时，提取科诺安图片地址
	sample_dict=[]
	find=False
	if len(exam_gd['clinic']):
		for item in exam_gd['clinic']:
			if item['bus_code']==sample_code.strip():
				sample_dict=item 
				find=True
	if find==False:
		print "The bus code:%s not found, please check the mongo url"%(sampleCode)
		sys.exit(1)
	exam_gd['imagepath']=sample_dict['imagepath'].strip()
	exam_gd['verify_cnv_data']=get_cnv_picture_result(bus_code=sample_code,status=exam_gd['imagepath'])
	exam_gd['verify_cnv_imgs']=trans_cnv_imgs(exam_gd['verify_cnv_data'])
	del exam_gd['verify_cnv_data']
	return exam_gd

def trans_vus_name(list):
	result=[]
	if len(list):
		for item in list:
			tmp=item
			if 'chr_position' in item.keys():
				tmp['chr_position']=item['chr_position']
			if 'mut_size' in item.keys():
				tmp['mut_size']=item['mut_size']
			if 'disease' in item.keys():
				tmp['Disease']=item['disease']
				del(tmp['disease'])
			if 'sample' in item.keys():
				tmp['Sample']=item['sample']
				del(tmp['sample'])
			if 'mut_asses' in item.keys():
				tmp['mut_asses']=item['mut_asses']
			if 'result' in item.keys():
				tmp['result']=item['result']
			if 'bus_code' in item.keys():
				tmp['bus_code']=item['bus_code']
			if 'positive_reason' in item.keys():
				tmp['positive_reason']=item['positive_reason']
			if 'variant_type' in item.keys():
				tmp['variant_type']=item['variant_type']
			if 'type' in item.keys():
				tmp['type']=item['type']
			if 'gene_list' in item.keys():
				tmp['gene_list']=item['gene_list']
			result.append(tmp)
	
	return result
	
def trans_note_name(list):
	result=[]
	if len(list):
		for item in list:
			tmp=item
			tmp['Sample']=item['sample']
			del(tmp['sample'])
			if 'gene' in item.keys():
				tmp['Gene']=item['gene']
				del(tmp['gene'])
			if 'position' in item.keys():
				tmp['Position']=item['position']
				del(tmp['position'])
			if 'NM' in item.keys():
				tmp['Transcript']=item['NM']
				del(tmp['NM'])
			if 'exon' in item.keys():
				tmp['Exon']=item['exon']
				del(tmp['exon'])
			if 'hgvs_c' in item.keys():
				tmp['hgvs.c']=item['hgvs_c']
				del(tmp['hgvs_c'])
			if 'hgvs_p' in item.keys():
				tmp['hgvs.p']=item['hgvs_p']
				del(tmp['hgvs_p'])
			if 'var_type' in item.keys():
				tmp['Type']=item['var_type']
				del(tmp['var_type'])
			if 'family_carry' in item.keys():
				tmp['Genotype']=item['family_carry']
				del(tmp['family_carry'])
				del(tmp['genotype'])
			if 'clinical_level' in item.keys():
				tmp['ACMGLevel']=item['clinical_level']
				del(tmp['clinical_level'])
			if 'disease' in item.keys():
				tmp['Disease']=item['disease']
				del(tmp['disease'])
			if 'inheritance' in item.keys():
				tmp['Inheritance']=item['inheritance']
				del(tmp['inheritance'])
			if 'maf' in item.keys():
				tmp['AlleleFrequency']=item['maf']
				del(tmp['maf'])
			if 'rs' in item.keys():
				tmp['dbSNP']=item['rs']
				del(tmp['rs'])
			result.append(tmp)
	
	return result
	
def trans_core_name(list):
	result=[]
	if len(list):
		for item in list:
			tmp=item
			tmp['Sample']=item['sample']
			del(tmp['sample'])
			tmp['Gene']=item['gene']
			del(tmp['gene'])
			tmp['Position']=item['position']
			del(tmp['position'])
			tmp['Transcript']=item['NM']
			del(tmp['NM'])
			tmp['Exon']=item['exon']
			del(tmp['exon'])
			tmp['hgvs.c']=item['hgvs_c']
			del(tmp['hgvs_c'])
			tmp['hgvs.p']=item['hgvs_p']
			del(tmp['hgvs_p'])
			tmp['Type']=item['var_type']
			del(tmp['var_type'])
			tmp['Genotype']=item['family_carry']
			del(tmp['family_carry'])
			tmp['ACMGLevel']=item['clinical_level']
			del(tmp['clinical_level'])
			tmp['Disease']=item['disease']
			del(tmp['disease'])
			tmp['Inheritance']=item['inheritance']
			del(tmp['inheritance'])
			tmp['AlleleFrequency']=item['maf']
			del(tmp['maf'])
			tmp['dbSNP']=item['rs']
			del(tmp['rs'])
			result.append(tmp)
	return result
def trans_verify_imgs(dict):
	index=1
	verify_imgs={}
	for item in dict:
		img_title=dict[item][1].split("/")[-1]
		verify_imgs[index]=[dict[item][0],img_title]
		index=index+1
	return verify_imgs
def trans_cnv_imgs(list):
	verify_cnv_img={}
	cnv_imgs={}
	imgs=[]
	for item in list:
		img_title=item.split("/")[-1]
		imgs.append(img_title)
	cnv_imgs=(',').join(imgs)
	return cnv_imgs

	
def  main():
	usage = "Usage: %prog -i input_data_file -b busCode -a ana_type -c check_result -o output "
	parser=OptionParser(usage)
	parser.add_option("-i", "--input", dest="input", action="store", help="txt file contains the result information")
	parser.add_option("-b","--buscode",dest="buscode",action="store",help="The sample Business Code")
	parser.add_option("-a","--ana_type",dest="ana_type",help="The analysis type,option: bulk or family")
	parser.add_option("-c","--check_result",dest="check_result",help="The sample check result, option:positive or negative")
	parser.add_option("-o","--output",dest="output",help="The output directory name")
	parser.add_option("-y","--hospital",dest="hospital",default="jiahui",help="which hospital need the report, if the hospital is xinhua hospital,the parameter needed!")
	(options,args)=parser.parse_args()
	global pdf_name
	output_pdf_name=""
	all_json_file=""
	if not options.input:  
		print "The input data file contains the result information must be specified!"
		sys.exit(1)
	if not options.buscode:
		print "The business code of the sample must be specified!"
		sys.exit(2)
	if not options.ana_type or (options.ana_type!="bulk" and options.ana_type!="family"):
		print "The sample analysis type must be specified, options: bulk or family!"
		sys.exit(4)
	if not options.check_result or (options.check_result!="positive" and options.check_result!="negative"):
		print "The sample check result must be specified, options: negative or positive!"
		sys.exit(5)	
	if not options.output:
		output_pdf_name=options.buscode+".检测报告.pdf" 
		all_json_file=absolute_json+options.buscode+".sup.json"###新增
	if options.output:
		if not os.path.exists(options.output):
			os.makedirs(options.output)
		output_pdf_name=options.buscode+".检测报告.pdf" 
		all_json_file=absolute_json+options.buscode+".sup.json"###新增
	###新增命令调用输出
	CMD_path=os.path.abspath(sys.argv[0])
	input_path_prefix='/share/production/Genetics/WES/wes_website/Upload/*/'
	sample_input_path=input_path_prefix+options.buscode+".report.conf"
	abstract_path=glob.glob(sample_input_path)
	CMD='python2 '+CMD_path+' -i '+ ''.join(abstract_path) +' -b '+ options.buscode
	if options.ana_type=="bulk" and options.check_result=="positive":
		#----1 读取输入模板，转换成xml格式，并获取结论、结果解读、一代验证图片、cnv图片、文献等信息
		create_common_xml.buildNewsXmlFile(options.input)
		import read_xml
		supplemnt = read_xml.supplement_dictionary()
		##---1.1 提取一代验证真正的图信息,去掉WES检测范围内的CNV图片
		#-------1.1 提取一代验证真正的图信息,以及WES检测范围内的CNV图片
		wes_verify_imgs,wescnv_imgs_names=get_verify_imgs_wescnv_imgs(supplemnt['verify_result'])
		#----------1.1.1实际的一代验证图片####
		supplemnt['verify_result']=wes_verify_imgs
		##---2 获取gd中补充位点：supplement_report和negative_appendix的数据##############
		exam_gd=get_gd_info(options.buscode,all_json_file)
		##---2.1 获取补充类型、补充原因#############################################
		###暂定补充类型和原因一致modify 20190327
		exam_gd['supplement_info']=supplemnt['supplement_info']
		
		###如果为加验、改做WES+再分析时，文献为空
		if exam_gd['supplement_type']=='add_sup_sites' or exam_gd['supplement_type']=='sup_ne_appendix':
			exam_gd['articles']=''
		else:
			exam_gd['articles']=supplemnt['articles']
		###阴转阳的位点解读是存在的，否则为空，且需要标注红色，加粗
		exam_gd['verify_result'] = trans_verify_imgs(wes_verify_imgs)
		exam_gd['verify_wescnv_imgs'] = wescnv_imgs_names
		exam_gd['cnv_seq']=supplemnt['cnv_seq']
		exam_gd['conclusion_summary']=supplemnt['conclusion_summary']		
		exam_gd['red']=supplemnt['red']#add20190214
		exam_gd['overstriking']=supplemnt['overstriking']#add20190214
		
		##---3 统一字段名称
		##增加补充报告的字段统一，add20190325同note
		if len(exam_gd['supplement_report']):
			exam_gd['supplement_sites']=trans_note_name(exam_gd['supplement_report'])	
		elif len(exam_gd['negative_appendix']):
			exam_gd['supplement_sites']=trans_ne_appendix_name(exam_gd['negative_appendix'])
		del exam_gd['supplement_report']
		del exam_gd['negative_appendix']
		
		##---4 将NA改为ND
		##增加补充报告NA转为ND，add20190325同note
		exam_gd['supplement_sites']=trans_NA_to_ND(exam_gd['supplement_sites'])
		##---5 去掉不需要的数据
		del exam_gd['result']
		del exam_gd['note']
		del exam_gd['vus']
		del exam_gd['clinic']
		##---6 封装json############
		data_string=json.dumps(exam_gd)
		f=open(all_json_file,"w")
		f.write(data_string)
		f.write('\n')
		##---7 要把json字符串post给IT的生成报告
		print "开始对接 : %s" % time.ctime()
		print "传递的JSON如下："
		if exam_gd['subnumber']=="":
			print "缺少子case编号！"
			sys.exit(1)
		''''url="http://eipapi.berrygenomics.com/berry/Jrest/wesAnlsReceiveService/addWesSupAnlsData"
		headers = {'Content-Type': 'application/json'}
		request = urllib2.Request(url=url, headers=headers, data=json.dumps(exam_gd))
		response = urllib2.urlopen(request).read()
		print "返回值：",response
		'''
		##对接CMD
		if '阳性' in exam_gd['check_result']:
			CMD=CMD+' -a bulk -c '+'positive'
		print '对接程序：',CMD
		print '\n'
		print "对接结束 : %s" % time.ctime()
		#supplement_bulk_postive_data(exam_dict, output_pdf_name,supplement_dict=new_supplemnt)
	elif options.ana_type=="bulk" and options.check_result=="negative":
		#----1 读取输入模板，转换成xml格式，并获取结论、结果解读、一代验证图片、cnv图片、文献等信息
		create_common_xml.buildNewsXmlFile(options.input)
		import read_xml
		supplemnt = read_xml.supplement_dictionary()
		##---1.1 提取一代验证真正的图信息,去掉WES检测范围内的CNV图片
		wes_verify_imgs,wescnv_imgs_names=get_verify_imgs_wescnv_imgs(supplemnt['verify_result'])
		#----------1.1.1实际的一代验证图片####
		supplemnt['verify_result']=wes_verify_imgs
		########################更新获取json的方式
		##---2 获取gd中补充位点：supplement_report和negative_appendix的数据##############
		exam_gd=get_gd_info(options.buscode,all_json_file)
		##---2.1 获取补充类型、补充原因#############################################
		exam_gd['supplement_info']=supplemnt['supplement_info']
		
		###转换一代验证
		exam_gd['verify_result'] = trans_verify_imgs(wes_verify_imgs)
		exam_gd['verify_wescnv_imgs'] = wescnv_imgs_names
		###如果为加验、改做WES+再分析时，文献为空
		if exam_gd['supplement_type']=='add_sup_sites' or exam_gd['supplement_type']=='sup_ne_appendix':
			exam_gd['articles']=''
		else:
			exam_gd['articles']=supplemnt['articles']
		###阴转阳的位点解读是存在的，否则为空，且需要标注红色，加粗
		exam_gd['cnv_seq']=supplemnt['cnv_seq']
		exam_gd['conclusion_summary']=supplemnt['conclusion_summary']		
		exam_gd['red']=supplemnt['red']#add20190214
		exam_gd['overstriking']=supplemnt['overstriking']#add20190214		
		##---3 统一字段名称
		if len(exam_gd['supplement_report']):
			exam_gd['supplement_sites']=trans_note_name(exam_gd['supplement_report'])	
		elif len(exam_gd['negative_appendix']):
			exam_gd['supplement_sites']=trans_ne_appendix_name(exam_gd['negative_appendix'])
		del exam_gd['supplement_report']
		del exam_gd['negative_appendix']
		##---4 将NA改为ND
		exam_gd['supplement_sites']=trans_NA_to_ND(exam_gd['supplement_sites'])
		##---5 去掉不需要的数据
		del exam_gd['result']
		del exam_gd['note']
		del exam_gd['vus']
		del exam_gd['clinic']
		##---6 封装json############
		data_string=json.dumps(exam_gd)
		f=open(all_json_file,"w")
		f.write(data_string)
		f.write('\n')
		if exam_gd['subnumber']=="":
			print "缺少子case编号！"
			sys.exit(1)
		##---7 要把json字符串post给IT的生成报告
		print "开始对接 : %s" % time.ctime()
		print "传递的JSON如下："
		'''url="http://eipapi.berrygenomics.com/berry/Jrest/wesAnlsReceiveService/addWesSupAnlsData"
		headers = {'Content-Type': 'application/json'}
		request = urllib2.Request(url=url, headers=headers, data=json.dumps(exam_gd))
		response = urllib2.urlopen(request).read()
		print "返回值：",response
		'''
		##对接CMD
		if '阴性' in exam_gd['check_result']:
			CMD=CMD+' -a bulk -c '+'negative'
		if '未知' in exam_gd['check_result']:
			CMD=CMD+' -a bulk -c '+'unknown'
		print '对接程序：',CMD
		print '\n'
		print "对接结束 : %s" % time.ctime()
	elif options.ana_type=="family" and options.check_result=="negative":
		#----1 读取输入模板，转换成xml格式，并获取结论、结果解读、一代验证图片、cnv图片、文献等信息
		create_family_xml.buildNewsXmlFile(options.input)
		import read_family_xml
		#exam_dict=get_sample_info(options.buscode, modify_result={'check_result': u"阴性", "ana_type": u"家系"},hospital=options.hospital)
		supplemnt=read_family_xml.supplement_dictionary()
		##---1.1 提取一代验证真正的图信息,去掉WES检测范围内的CNV图片
		#-------1.1 提取一代验证真正的图信息,以及WES检测范围内的CNV图片
		wes_verify_imgs,wescnv_imgs_names=get_verify_imgs_wescnv_imgs(supplemnt['verify_result'])
		#----------1.1.1实际的一代验证图片####
		supplemnt['verify_result']=wes_verify_imgs
		########################更新获取json的方式
		##---2 获取gd 位点数据及qc，extend路径的数据##############
		exam_gd=get_gd_info(options.buscode,all_json_file)
		##---2.1 获取补充类型、补充原因#############################################
		exam_gd['supplement_info']=supplemnt['supplement_info']
		###转换一代验证
		exam_gd['verify_result'] = trans_verify_imgs(wes_verify_imgs)
		exam_gd['verify_wescnv_imgs'] = wescnv_imgs_names
		###如果为加验、改做WES+再分析时，文献为空
		if exam_gd['supplement_type']=='add_sup_sites' or exam_gd['supplement_type']=='sup_ne_appendix':
			exam_gd['articles']=''
		else:
			exam_gd['articles']=supplemnt['articles']
		###阴转阳的位点解读是存在的，否则为空，且需要标注红色，加粗
		exam_gd['cnv_seq']=supplemnt['cnv_seq']
		exam_gd['conclusion_summary']=supplemnt['conclusion_summary']	
		exam_gd['red']=supplemnt['red']#add20190214
		exam_gd['overstriking']=supplemnt['overstriking']#add20190214
		
		##---3 统一字段名称
		##增加补充报告的字段统一，add20190325同note
		if len(exam_gd['supplement_report']):
			exam_gd['supplement_sites']=trans_note_name(exam_gd['supplement_report'])	
		elif len(exam_gd['negative_appendix']):
			exam_gd['supplement_sites']=trans_ne_appendix_name(exam_gd['negative_appendix'])
			
		del exam_gd['supplement_report']
		del exam_gd['negative_appendix']
		##---4 将NA改为ND
		exam_gd['supplement_sites']=trans_NA_to_ND(exam_gd['supplement_sites'])
		##---5 去掉不需要的数据
		del exam_gd['result']
		del exam_gd['note']
		del exam_gd['vus']
		del exam_gd['clinic']
		
		##---7 封装json############
		data_string=json.dumps(exam_gd)
		f=open(all_json_file,"w")
		f.write(data_string)
		f.write('\n')
		##---8 要把json字符串post给IT的生成报告
		print "开始对接 : %s" % time.ctime()
		print "传递的JSON如下："
		print "\n"
		if exam_gd['subnumber']=="":
			print "缺少子case编号！"
		'''	sys.exit(1)
		url="http://eipapi.berrygenomics.com/berry/Jrest/wesAnlsReceiveService/addWesSupAnlsData"
		headers = {'Content-Type': 'application/json'}
		request = urllib2.Request(url=url, headers=headers, data=json.dumps(exam_gd))
		response = urllib2.urlopen(request).read()
		print "返回值：",response
		'''
		####CMD#######
		if '阴性' in exam_gd['check_result']:
			CMD=CMD+' -a family -c '+'negative'
		if '未知' in exam_gd['check_result']:
			CMD=CMD+' -a family -c '+'unknown'
		print '对接程序：',CMD
		print '\n'
		print "对接结束 : %s" % time.ctime()
	elif options.ana_type=="family" and options.check_result=="positive":
		#----1 读取输入模板，转换成xml格式，并获取结论、结果解读、一代验证图片、cnv图片、文献等信息
		create_family_xml.buildNewsXmlFile(options.input)
		import read_family_xml  
		#exam_dict=get_sample_info(options.buscode,modify_result={'check_result': u"阳性", "ana_type": u"家系"},hospital=options.hospital)
		supplemnt=read_family_xml.supplement_dictionary()
		##---1.1 提取一代验证真正的图信息,去掉WES检测范围内的CNV图片
		wes_verify_imgs,wescnv_imgs_names=get_verify_imgs_wescnv_imgs(supplemnt['verify_result'])
		#----------1.1.1实际的一代验证图片####
		supplemnt['verify_result']=wes_verify_imgs
		##---2 获取gd 位点数据及qc，extend路径的数据##############
		exam_gd=get_gd_info(options.buscode,all_json_file)
		exam_gd['verify_wescnv_imgs']=wescnv_imgs_names
		####2.1 获取家系其他成员的CNV信息,合并
		exam_gd['vus']=get_other_family_members_cnv(options.buscode,exam_gd['bus_code_combine'],exam_gd['vus'])
		##---2.1 获取补充类型、补充原因#############################################
		exam_gd['supplement_info']=supplemnt['supplement_info']
		###转换一代验证
		exam_gd['verify_result'] = trans_verify_imgs(wes_verify_imgs)
		exam_gd['verify_wescnv_imgs'] = wescnv_imgs_names
		###如果为加验、改做WES+再分析时，文献为空
		if exam_gd['supplement_type']=='add_sup_sites' or exam_gd['supplement_type']=='sup_ne_appendix':
			exam_gd['articles']=''
		else:
			exam_gd['articles']=supplemnt['articles']
		###阴转阳的位点解读是存在的，否则为空，且需要标注红色，加粗
		exam_gd['cnv_seq']=supplemnt['cnv_seq']
		exam_gd['conclusion_summary']=supplemnt['conclusion_summary']	
		exam_gd['red']=supplemnt['red']#add20190214
		exam_gd['overstriking']=supplemnt['overstriking']#add20190214
		
		##---3 统一字段名称
		##增加补充报告的字段统一，add20190325同note
		if len(exam_gd['supplement_report']):
			exam_gd['supplement_sites']=trans_note_name(exam_gd['supplement_report'])	
		elif len(exam_gd['negative_appendix']):
			exam_gd['supplement_sites']=trans_ne_appendix_name(exam_gd['negative_appendix'])
		del exam_gd['supplement_report']
		del exam_gd['negative_appendix']
		##---4 将NA改为ND
		exam_gd['supplement_sites']=trans_NA_to_ND(exam_gd['supplement_sites'])
		##---5 去掉不需要的数据
		del exam_gd['result']
		del exam_gd['note']
		del exam_gd['vus']
		del exam_gd['clinic']
		##---7 封装json############
		data_string=json.dumps(exam_gd)
		f=open(all_json_file,"w")
		f.write(data_string)
		f.write('\n')
		##---8 要把json字符串post给IT的生成报告
		print "开始对接 : %s" % time.ctime()
		print "传递的JSON如下："
		if exam_gd['subnumber']=="":
			print "缺少子case编号！"
			sys.exit(1)
		'''url="http://eipapi.berrygenomics.com/berry/Jrest/wesAnlsReceiveService/addWesSupAnlsData"
		headers = {'Content-Type': 'application/json'}
		request = urllib2.Request(url=url, headers=headers, data=json.dumps(exam_gd))
		response = urllib2.urlopen(request).read()
		print "返回值：",response
		'''
		####CMD#######
		if '阳性' in exam_gd['check_result']:
			CMD=CMD+' -a family -c '+'positive'
		print '对接程序', CMD
		print '\n'
		print "对接结束 : %s" % time.ctime()
										
if	__name__=="__main__":
	main()   

