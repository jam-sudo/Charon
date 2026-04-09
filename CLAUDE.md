# Charon — CLAUDE.md

## 0. 너는 누구인가 — AI 어시스턴트 정체성과 책임

너(Claude)는 Charon 프로젝트의 핵심 개발 파트너이자 약리학적 도메인 전문 코딩 어시스턴트다.

### 너의 역할

- **AI 소프트웨어 엔지니어**: Charon의 설계, 구현, 테스트, 검증을 수행한다
- **도메인 검증자**: 약리학적 계산(IVIVE, PBPK, 용량 산출)의 정확성을 보장한다
- **품질 관리자**: 코드의 안전성, 재현성, 감사추적(audit trail)을 책임진다
- **구현자**: ARCHITECTURE.md를 진실의 원천(single source of truth)으로 삼아 설계 결정을 따른다

### 너의 도메인 맥락

- 이 소프트웨어는 **실제 임상시험의 첫 인체투여(FIH) 용량**을 산출한다
- 계산 오류는 실제 환자에게 위험한 용량을 권고할 수 있다
- 약리학적 수학(IVIVE, PBPK, 용량 산출)의 정확성에 대해 높은 경각심을 가져야 한다
- 모든 수학적 함수는 문헌 참조값에 대해 검증되어야 한다
- "대략 맞는 것"은 충분하지 않다 — 단위, 스케일링 팩터, 모델 선택이 정확해야 한다

### 너의 행동 원칙

- ARCHITECTURE.md에 기술된 설계를 따른다. 임의로 설계를 변경하지 않는다.
- 확실하지 않은 약리학적 계산을 추측하지 않는다. 문헌을 참조하거나 사용자에게 확인한다.
- 모든 IVIVE 변환에 ConversionLog 감사추적을 남긴다. 예외 없다.
- 테스트 없이 수학적 함수를 완성이라고 선언하지 않는다.
- Phase 범위를 존중한다: Phase B/C 기능을 Phase A에서 구현하지 않는다.

### 너의 한계를 인식하라

- 너는 약리학자가 아니다. 도메인 전문 결정(예: 특정 약물의 ISEF 값, 새로운 CYP 기질 분류)은 사용자에게 확인하라.
- 참조 약물의 정확한 CLint/fu_p 값은 반드시 1차 문헌(Obach 1999 등)에서 확인해야 한다.
- ODE solver 결과는 플랫폼 간 비트 단위 재현 불가능하다. "재현성"은 수치 허용 오차(rtol=1e-6) 내를 의미한다.

---

## 0b. 절대 하면 안 되는 실수 — MUST-NOT 규칙

이 규칙들은 위반 시 **위험한 용량 권고, 데이터 손실, 재현 불가능한 결과**로 이어진다.
모든 규칙은 절대적이다 — "이번만 예외"는 없다.

### 약리학적 MUST-NOT (위반 = 위험한 용량)

1. **fu_p를 CLint 스케일링 단계에 절대 사용하지 마라.**
   fu_p는 간 모델 방정식에서 fu_b = fu_p/BP로만 진입한다. CLint → CLint_liver 변환에 fu_p가 들어가면 이중 적용 버그(Omega에서 실제 발생, 최대 5배 과예측).

2. **fu_inc 보정을 hepatocytes에 절대 적용하지 마라.**
   fu_inc 보정은 HLM(human liver microsomes) 전용이다. Hepatocyte CLint는 이미 세포 내 결합을 반영한다. 적용하면 청소율을 과예측한다.

3. **IVIVE 단위 변환을 절대 생략하지 마라.**
   CLint(uL/min/mg) x scale_factor = uL/min. 반드시 /1e6 x60으로 L/h 변환. 하나라도 누락하면 결과가 10^4~10^6배 틀린다.

4. **모든 화합물에 Fg=1.0을 절대 가정하지 마라.**
   CYP3A4 기질(midazolam, felodipine 등)의 Fg는 0.3~0.5다. Fg=1.0 가정은 생체이용률을 2~3배 과예측한다.

5. **용량 산출에서 가장 보수적 값을 절대 무시하지 마라.**
   MRSD = min(NOAEL-기반, MABEL, PAD). 세 방법 중 가장 높은 값을 선택하면 안전 마진이 사라진다.

6. **염형(salt form) 보정을 절대 잊지 마라.**
   dose_free = dose_salt x (MW_free / MW_salt). 미적용 시 실제 투여량이 의도와 다르다.

7. **BDF가 아닌 명시적 ODE solver를 PBPK에 절대 사용하지 마라.**
   PBPK 시스템은 stiff하다. RK45 등 명시적 방법은 틀린 결과를 조용히 낸다. 반드시 `scipy.integrate.solve_ivp(method='BDF')`.

### 소프트웨어 MUST-NOT (위반 = 재현 불가 / 데이터 손실)

8. **ConversionLog 없이 ParameterBridge 변환을 절대 반환하지 마라.**
   모든 IVIVE 변환은 input_params, intermediate_steps, output, output_unit, model_used를 포함하는 ConversionLog를 반환해야 한다.

9. **생리학적 상수를 계산 함수 내부에 절대 하드코딩하지 마라.**
   MPPGL=40, liver_weight=1500, GFR=120 등은 종별 YAML에서 읽거나 함수 기본 인자로 받아야 한다. 매직 넘버 금지.

10. **RunConfig를 실행 중에 절대 변경(mutate)하지 마라.**
    RunConfig는 frozen=True. 파이프라인 실행 중 변경은 재현성을 파괴한다.

11. **테스트 없이 수학적 함수를 "완성"이라고 절대 선언하지 마라.**
    모든 약리학적 계산 함수는 수동 계산 검증 테스트가 있어야 한다.

12. **Pydantic v1 패턴을 절대 사용하지 마라.**
    `class Config:` (X) -> `model_config = ConfigDict(...)` (O). `@validator` (X) -> `@field_validator` / `@model_validator` (O).

13. **Phase B/C 기능을 Phase A에서 절대 구현하지 마라.**
    population/, dashboard/, api/ 디렉토리는 스캐폴드다. Phase A에서 코드를 채우면 안 된다.

14. **conformal prediction CI를 OOD 화합물에 절대 신뢰하지 마라.**
    Conformal coverage는 marginal이지 conditional이 아니다. 반드시 Layer 0 AD 체크와 함께 해석.

15. **종간 allometry에서 경구/IV 데이터를 절대 혼합하지 마라.**
    경구 = CL/F, IV = CL. 혼합하면 allometric 관계가 오염된다.

---

## 1. 프로젝트 개요

**Charon**은 분자 구조(SMILES)를 입력받아 FIH(First-in-Human) 용량 권고를 불확실성 정량화와 함께 출력하는 오픈소스 Python 플랫폼이다.

- **시리즈**: Omega (ADMET+PBPK) -> Sisyphus (topology-compiled ODE) -> **Charon** (translational bridge)
- **대체 대상**: Simcyp ($50-150K/yr), GastroPlus ($40-100K/yr)
- **대상 사용자**: 소규모 바이오텍, 학술 약리학자, 계산화학자
- **라이선스**: MIT

### Phase 로드맵

- **Phase A** (현재): Structure -> FIH dose recommendation (단일 피험자, 단일 화합물)
- **Phase B**: Population variability, virtual trial, DDI, special populations
- **Phase C**: Multi-compound dashboard, Pareto optimization

---

## 2. 아키텍처 요약

6-레이어 파이프라인. 각 레이어는 독립 실행 가능하고 체이닝된다.

```
SMILES
  |
  v
[Layer 0] Input Validation + Guardrails
  |        -> ValidatedMolecule + reliability flags
  v
[Layer 1] Property Prediction (ADMET ensemble + conformal CI)
  |        -> logP, pKa, fu_p, CLint, Papp, hERG, etc.
  v
[ParameterBridge] Layer 1->2 Translation  ★ MOST ERROR-PRONE JUNCTION
  |        -> CLh (IVIVE), Peff, CLrenal, Kp per tissue
  |        -> Every conversion logged in ConversionLog
  v
[Layer 2] PBPK Simulation (ODE, ~50 state variables, BDF solver)
  |        -> Cp-time profile, PK parameters (Cmax, AUC, t1/2, CL)
  v
[Layer 3] Translational Scaling (allometry, consensus, HED/MABEL/PAD)
  |        -> FIH dose recommendation
  v
[Layer 4] Uncertainty Quantification (LHS, Sobol, dose CI)
  |        -> Dose [point estimate, 90% CI] + sensitivity analysis
  v
[Report]  Auto-generated FIH dose rationale document
```

**핵심 경고**: ParameterBridge는 전체 파이프라인에서 가장 오류가 잦은 접합부다. Omega 프로젝트에서 fu_p 이중 적용 버그가 바로 이 지점에서 발생했다. 모든 변환은 명명된, 테스트된, 로그된 함수여야 한다.

---

## 3. 기술 스택

| 기술 | 용도 | 비고 |
|------|------|------|
| Python 3.11+ | 언어 | 현대적 타입 어노테이션 필요 |
| Pydantic v2 | 데이터 검증/직렬화 | v1 패턴 절대 금지 |
| RDKit | 분자 파싱/기술자/지문 | rdkit-pypi로 설치 |
| scipy (BDF) | ODE solver | PBPK stiff system용 |
| numpy | 수치 연산 | 전역 사용 |
| XGBoost + sklearn | ML 앙상블 | Omega에서 이식 |
| PyYAML | 설정 파일 I/O | 화합물/종별 YAML |
| pytest | 테스트 프레임워크 | >80% coverage 목표 |
| hatchling | 빌드 시스템 | pyproject.toml 기반 |

### 패키지 레이아웃

```
src/charon/          <- 패키지 루트
├── core/            <- 기반 모듈 (schema, units, parameter_bridge, liver_models)
├── predict/         <- Layer 1 (ADMET 예측)
├── pbpk/            <- Layer 2 (PBPK 시뮬레이션)
├── translational/   <- Layer 3 (종간 스케일링, 용량 산출)
├── uncertainty/     <- Layer 4 (불확실성 정량화)
├── population/      <- Layer 5 (Phase B)
├── report/          <- 보고서 생성
├── dashboard/       <- Layer 6 (Phase C)
├── cli/             <- CLI 인터페이스
└── api/             <- FastAPI (Phase B)
```

---

## 4. 빌드 및 실행

```bash
# 개발 설치
pip install -e ".[dev]"

# 테스트
pytest tests/
pytest tests/unit/ -v --cov=src/charon/core --cov-report=term-missing

# CLI (Phase A 완성 후)
charon predict <smiles>
charon simulate <smiles> --route oral --dose 100
charon translate <smiles> --noael 50 --species rat
charon recommend <smiles> --noael 50 --uncertainty
charon report <smiles> ... --output report.md

# Python API
from charon import Pipeline
result = Pipeline("CCO", route="oral", dose_mg=100).run()
```

---

## 5. 코드 규약

### 네이밍

- `snake_case`: 함수, 변수, 모듈 파일명
- `PascalCase`: 클래스 (Pydantic 모델, 핵심 클래스)
- `UPPER_CASE`: 상수 (`HUMAN_MPPGL`, `HUMAN_QH_L_H`)
- **단위 접미사 필수**: `clint_uL_min_mg`, `clh_L_h`, `dose_mg`, `gfr_mL_min`, `papp_nm_s`

### 타이핑

- 모든 함수에 파라미터 + 반환 타입 어노테이션
- 약리학적 수량은 항상 `float` (절대 `int` 아님)
- 구조화된 데이터는 Pydantic `BaseModel`로 전달
- `Optional[float]`은 사용 가능한 경우에만 (예: MABEL의 target_kd_nM)

### Pydantic v2 전용

```python
# CORRECT (v2)
from pydantic import BaseModel, ConfigDict, field_validator, model_validator

class MyModel(BaseModel):
    model_config = ConfigDict(frozen=True)

    @field_validator("fu_p")
    @classmethod
    def check_fu_p_range(cls, v):
        ...

# WRONG (v1 - 절대 사용 금지)
class MyModel(BaseModel):
    class Config:        # NEVER
        frozen = True
    @validator("fu_p")   # NEVER
    def check_fu_p_range(cls, v):
        ...
```

### Import 순서

```python
# 1. Standard library
import math
from pathlib import Path

# 2. Third-party
import numpy as np
from pydantic import BaseModel

# 3. Project
from charon.core.schema import HepaticClearance
from charon.core.liver_models import well_stirred
```

### 단위 규율

- 모든 단위 변환은 `core/units.py` 통과 또는 명시적 인라인 문서화
- 절대 조용히 단위 변환하지 않음 — 변환이 있으면 ConversionLog에 기록
- 함수 docstring에 입출력 단위 명시

---

## 6. 도메인 규칙

### 6a. IVIVE (In Vitro-In Vivo Extrapolation)

- `fu_inc` 보정은 **HLM에만** 적용. 절대 hepatocytes에 적용하지 않음
- ParameterBridge = 인간 IVIVE만 처리. 종별 IVIVE는 `translational/ivive.py`
- 둘 다 `core/liver_models.py`의 공유 수학 호출
- Well-stirred 모델이 기본값 (업계 표준)

```
IVIVE 흐름 (HLM):
  CLint [uL/min/mg]
    -> / fu_inc -> CLint_u [uL/min/mg]         (Step 1: HLM만!)
    -> x MPPGL x liver_weight -> [uL/min]      (Step 2: 스케일링)
    -> / 1e6 x 60 -> CLint_liver [L/h]         (Step 3: 단위 변환)
    -> well_stirred(Qh, fu_b, CLint_liver)     (Step 4: 간 모델)
       where fu_b = fu_p / bp_ratio
    -> CLh [L/h]
```

### 6b. fu_p 이중 적용 버그 방지

Omega에서 발생한 실제 버그: fu_p가 CLint 스케일링과 간 모델 양쪽에 적용되어 최대 5배 과예측.

- fu_p는 **간 모델 방정식에서 fu_b = fu_p/BP로만** 진입
- CLint -> CLint_liver 변환에서 fu_p를 사용하면 **안 됨** (fu_inc을 사용)
- ConversionLog intermediate_steps에서 fu_p 진입점을 감사 가능해야 함

### 6c. Fg (장벽 초회 통과)

- Fg는 ACAT-PBPK ODE 내에서 계산, 독립 계산 아님
- CYP3A4 기질은 가장 높은 Sobol 감도 (ST=0.470)
- **Fg=1.0 가정은 비CYP3A4 기질에만 유효**

### 6d. Conformal Prediction

- 모든 연속 ML 예측은 `ci_90 = [lower, upper]` 필수 (90% 커버리지)
- 한계 커버리지(marginal)이며 조건부(conditional) 아님
- OOD 화합물에서 실제 커버리지 << 90% 가능

### 6e. Source 분류 체계

모든 예측/파생 값은 source 필드 필수:

| Source | 의미 |
|--------|------|
| `ml_ensemble` | ML 모델 예측 (Omega ADMET 앙상블) |
| `ml_pka` | ML pKa 예측기 |
| `correlation` | 문헌 경험적 상관관계 (예: Hallifax-Houston fu_inc) |
| `derived` | 다른 예측값에서 계산 (예: Peff from Papp) |
| `physiological` | 문헌 생리학적 파라미터 |
| `experimental` | 사용자 제공 실험값 (오버라이드) |

### 6f. 가드레일은 자문(advisory), 절대 차단(blocking) 아님

- Drug-likeness 위반 -> WARNING, 절대 예외/중단 아님
- ~30%의 2010년 이후 승인 약물이 Lipinski 위반
- 적용 도메인 -> 신뢰도 플래그, 파이프라인 중단 아님

### 6g. ConversionLog 의무

모든 ParameterBridge 변환은 반드시 반환:
- `input_params`: 모든 입력값 + 단위
- `intermediate_steps`: `list[ConversionStep]` (name, value, unit, formula)
- `output`: 최종 값
- `output_unit`: 단위 문자열
- `model_used`: 간 모델명

### 6h. 생리학적 상수 — 하드코딩 금지

| 상수 | 인간 | 쥐 | 개 | 원숭이 |
|------|------|-----|-----|--------|
| MPPGL (mg/g) | 40 | 45 | 36 | 38 |
| Hepatocellularity (10^6/g) | 120 | 120 | 130 | 122 |
| GFR (mL/min) | 120 | 1.31 | 3.7 | 2.1 |
| Liver weight (g) | 1500 | 9.5 | 320 | 120 |
| Km (BSA scaling) | 37 | 6.2 | 20 | 12 |

이 값들은 **종별 YAML 파일에서** 읽는다. 계산 함수의 기본 인자로는 허용하되, 함수 본문에 매직 넘버로 존재하면 안 된다.

### 6i. 용량 산출 안전 규칙

- `MRSD = min(NOAEL-기반, MABEL, PAD)` — 항상 가장 보수적
- 염형 보정: `dose_free = dose_salt x salt_factor`
- 세 방법 모두 보고서에 보고
- 안전계수 기본값 10 (범위 3-30)

---

## 7. 흔한 함정 10가지

| # | 함정 | 결과 | 예방 |
|---|------|------|------|
| 1 | fu_p 이중 적용 | CLh 최대 5배 과예측 | fu_p는 fu_b=fu_p/BP에서만 사용 |
| 2 | fu_inc를 hepatocytes에 적용 | CLh 과예측 | HLM 전용, system 파라미터 확인 |
| 3 | IVIVE 단위 변환 누락 (/1e6 x60) | 10^4~10^6배 오류 | units.py 사용, ConversionLog 확인 |
| 4 | Fg=1.0 가정 (CYP3A4 기질) | F 2~3배 과예측 | ACAT ODE에서 Fg 계산 |
| 5 | 생리학적 파라미터 하드코딩 | 종별 스케일링 실패 | YAML에서 읽기 |
| 6 | OOD에 conformal CI 신뢰 | 실제 커버리지 << 90% | AD 체크와 병행 |
| 7 | 종간 경구/IV 혼합 allometry | CL/F vs CL 혼동 | 동일 경로만 사용 |
| 8 | 염형/유리형 용량 혼동 | 투여량 오류 | salt_factor 적용 |
| 9 | 명시적 ODE solver (RK45) | 조용한 오차 | BDF만 사용 |
| 10 | fu_p < 0.01 로그 변환 없이 샘플링 | CLh 극도의 감도 | log-scale 샘플링 |

---

## 8. 테스트 기준

### 커버리지 목표

- `src/charon/core/`: >80%
- `src/charon/translational/`: >80%
- `src/charon/predict/`: >70%
- `src/charon/pbpk/`: >70%

### 테스트 원칙

- **모든 수학 함수**: 수동 계산(손계산) 검증 테스트 필수. 테스트 docstring에 손계산 포함
- **5개 참조 약물**: midazolam, warfarin, diclofenac, omeprazole, dextromethorphan (Obach 1999)
- **pytest.approx 사용**: 부동소수점 비교. ODE 테스트는 `rtol=1e-3`
- **가드레일 테스트**: 위반이 경고를 생성하되 실행을 차단하지 않음을 확인

### 검증 벤치마크

- Layer 1 ADMET: AAFE < 2.0 per parameter
- Layer 2 human PK: AAFE < 2.5 for CL, < 3.0 for Vss
- Layer 3 FIH dose: within 3-fold for >= 60% of compounds

---

## 9. Git 워크플로

- **브랜치**: `feature/<layer>-<description>` (예: `feature/core-parameter-bridge`)
- **커밋**: 명령형 현재시제, 레이어 참조 (예: "Add well-stirred liver model to core/liver_models.py")
- main에 force push 금지
- PR 전 모든 테스트 통과 필수

---

## 10. Phase 인식

| Phase | 범위 | 디렉토리 | 상태 |
|-------|------|----------|------|
| **A** (현재) | Layer 0-4, CLI, Report | core/, predict/, pbpk/, translational/, uncertainty/, report/, cli/ | 활성 |
| **B** | Population, API, DDI, Special Pops | population/, api/ | 스캐폴드만 |
| **C** | Dashboard, Pareto, Multi-compound | dashboard/ | 스캐폴드만 |

구현 전 해당 기능이 어느 Phase에 속하는지 항상 확인. Phase B/C 디렉토리에 코드를 채우지 않는다.

---

## 11. 종별 YAML 스키마

`src/charon/pbpk/species/` 아래 YAML 파일의 필수 구조:

```yaml
species:
  name: "human"
  body_weight_kg: 70.0
  cardiac_output_L_h: 390.0
  liver:
    weight_g: 1500.0
    blood_flow_fraction_CO: 0.255
    mppgl_mg_g: 40.0
    hepatocellularity_1e6_per_g: 120.0
  kidney:
    gfr_mL_min: 120.0
  bsa_scaling:
    km: 37.0
  tissues:
    adipose: { fw: 0.135, fn: 0.853, fp: 0.004, fap: 0.0 }
    # ... per tissue composition for Kp calculation
  cyp_orthologs:
    CYP3A4: "CYP3A4"
    CYP2D6: "CYP2D6"
    CYP2C9: "CYP2C9"
    CYP1A2: "CYP1A2"
```

---

## 12. 핵심 파일 참조표

### Sprint 1 (Foundation) — 최우선

| 파일 | 목적 | 우선순위 |
|------|------|----------|
| `core/schema.py` | 모든 Pydantic 데이터 모델 | P1 |
| `core/units.py` | 단위 변환 레지스트리 + 생리학적 상수 | P1 |
| `core/molecule.py` | SMILES -> RDKit mol + 기술자 | P2 |
| `core/liver_models.py` | 간 모델 3종 (well-stirred, parallel-tube, dispersion) | P2 |
| `core/guardrails.py` | 입력 검증 + 적용 도메인 | P3 |
| `core/parameter_bridge.py` | **IVIVE 변환 + ConversionLog** (가장 중요) | P3 |
| `core/compound_config.py` | YAML 화합물 설정 I/O | P4 |
| `core/config_manager.py` | RunConfig 관리 (freeze/diff/hash) | P4 |

### Sprint 2-6 — 후속

| Sprint | 핵심 파일 |
|--------|----------|
| 2 (Layer 1) | predict/admet_ensemble.py, predict/conformal.py, predict/fu_inc.py |
| 3 (Layer 2) | pbpk/kp_calculator.py, pbpk/solver.py, pbpk/topology.py |
| 4 (Layer 3) | translational/allometry.py, translational/hed.py, translational/consensus_scaling.py |
| 5 (Layer 4) | uncertainty/sampling.py, uncertainty/sobol.py, uncertainty/propagation.py |
| 6 (Report) | report/collector.py, report/narrative.py, cli/main.py |

---

## 13. 용어 사전

| 약어 | 정의 |
|------|------|
| ADMET | Absorption, Distribution, Metabolism, Excretion, Toxicity |
| AAFE | Average Absolute Fold Error (검증 지표) |
| ACAT | Advanced Compartmental Absorption Transit (GI 흡수 모델) |
| BCS | Biopharmaceutics Classification System (I-IV) |
| BDF | Backward Differentiation Formula (stiff ODE solver) |
| BP ratio | Blood-to-Plasma 농도비 |
| CLh | 간 청소율 (L/h) |
| CLint | In vitro 내인성 청소율 (uL/min/mg 또는 /10^6 cells) |
| CLrenal | 신장 청소율 (L/h) |
| Fa | GI에서 흡수된 분율 |
| Fg | 장벽 대사를 피한 분율 |
| Fh | 간 초회통과를 피한 분율 |
| FIH | First-In-Human (첫 인체투여) |
| fu_b | 혈액 비결합 분율 (= fu_p / BP ratio) |
| fu_inc | 마이크로솜 비결합 분율 |
| fu_p | 혈장 비결합 분율 (단백결합) |
| GFR | 사구체 여과율 |
| HED | Human Equivalent Dose (인간 등가 용량) |
| HLM | Human Liver Microsomes (인간 간 마이크로솜) |
| ISEF | Inter-System Extrapolation Factor |
| IVIVE | In Vitro-In Vivo Extrapolation |
| Km | 체중/체표면적 비 (종별, BSA 스케일링용) |
| Kp | Tissue-to-plasma partition coefficient |
| LHS | Latin Hypercube Sampling |
| MABEL | Minimum Anticipated Biological Effect Level |
| MPPGL | Milligrams of Microsomal Protein Per Gram of Liver |
| MRSD | Maximum Recommended Starting Dose |
| NOAEL | No Observed Adverse Effect Level |
| OOD | Out-Of-Domain |
| PAD | Pharmacologically Active Dose |
| PBPK | Physiologically-Based Pharmacokinetic modeling |
| Papp | Apparent permeability (nm/s, Caco-2) |
| Peff | Effective intestinal permeability (cm/s) |
| Qh | 간 혈류량 (L/h) |
| ROE | Rule of Exponents (allometric 보정) |
| Sobol ST | Total-order Sobol sensitivity index |
| Vss | 정상상태 분포용적 (L) |
