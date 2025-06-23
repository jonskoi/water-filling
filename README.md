# 💧 Water-Filling Algorithm in OFDM Systems

OFDM 시스템에서는 여러 서브캐리어(subcarrier)를 사용해 동시에 데이터를 전송합니다.

하지만 무선 통신에서 **채널의 품질은 주파수마다 다르며**, 이는 곧 **서브캐리어마다 SNR이 다르다**는 것을 의미합니다.

> 📌 *문제*: **총 전력이 제한되어 있을 때**,  
> **어떻게 전력을 분배해야 전체 데이터 전송률(Capacity)을 극대화할 수 있을까요?**

---

## 🧠 왜 균등 분배는 비효율적일까?

서브캐리어마다 **채널 이득(Channel gain)** 이 다릅니다.

따라서 단순히 전력을 `1/N`로 똑같이 나누면 좋은 채널에도, 나쁜 채널에도 동일한 전력이 들어가 **비효율적**입니다.

이럴 때 사용하는 것이 바로 **워터필링(Water-Filling) 알고리즘**입니다.

---

## 💧 Water-Filling 이란?

“물을 붓는다”는 비유에서 이름이 유래되었습니다.

- 여러 개의 **높이가 다른 컵** = 서브캐리어들의 채널 역이득 (1/gain)
- 제한된 양의 **물** = 총 전력
- 각 컵의 높이에 따라 **물을 채우는 수준** = 전력 분배

> 낮은 컵(좋은 채널)부터 먼저 채우고, 나머지에 따라 높은 컵에도 채운다는 개념입니다.

---

## 🧮 최적화 문제 정의

총 N개의 서브캐리어에 대해, 서브캐리어 `k`의 채널 이득이 `gₖ`, 전력 할당량이 `εₖ`일 때 다음 최적화 문제를 풉니다:

**(최대화 목표)**

![Optimization Objective](https://github.com/user-attachments/assets/e5bef987-742d-4ccb-a5b2-dbb520c55866)

**제약 조건**

- $\sum_k \varepsilon_k \leq E$
- $\varepsilon_k \geq 0$

> 🔄 이 문제는 Concave Maximization → Convex Minimization으로 변환하여 풉니다.

![Convex Form](https://github.com/user-attachments/assets/2d75afc5-5cb3-4df3-ac90-f99a5ead7874)

---

## 🔧 KKT 조건을 통한 해 도출

Lagrangian 함수는 다음과 같이 정의됩니다:

![Lagrangian](https://github.com/user-attachments/assets/fc87aef8-1e89-435d-89bd-44cb2b10ec5c)

Lagrangian의 Gradient를 구하고:

![Gradient](https://github.com/user-attachments/assets/17a3b48f-746c-4052-b27d-b6c6b4d51f59)

**Complementary Slackness** 조건과 결합하여, 다음 식을 만족하는 값을 찾습니다:

![Complementary Slackness](https://github.com/user-attachments/assets/9e47de5f-c971-4974-b12e-c34093f9523b)

결국 각 서브캐리어의 최적 전력 할당은 다음과 같습니다:

![Final Solution](https://github.com/user-attachments/assets/f9709ca3-1e45-43d5-b62b-d1673ecc7708)

워터레벨 $\mu$는 다음과 같이 정의될 수 있습니다:

![Waterlevel](https://github.com/user-attachments/assets/8654743c-b35c-4f5c-ad5e-7016ab4284f3)

---

## 💻 MATLAB 시뮬레이션 코드

### ✅ Energy Allocation with Water-Filling

무작위 채널 이득에 대해 1000번 시뮬레이션하여 전력 분배 결과와 평균 전송률을 계산합니다.

```matlab
% 채널 이득 g와 노이즈 n에 대해 waterfilling 적용
% 전력 분배 결과 시각화 및 용량 비교
